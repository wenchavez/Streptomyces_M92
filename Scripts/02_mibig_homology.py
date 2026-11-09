from pathlib import Path
import pandas
import csv
import json
import subprocess
import re
import sys


base = Path(__file__).resolve().parent.parent

resultados = base/ "results"
salida = resultados/"MIBiG Homology"
sal_asignacion = salida/"structural allocation.csv"
sal_estructuras = salida/"putative structures.csv"


mibig_json = base/"Datos MIBiG"/"mibig_json_4.0"
mibig_gbk = base/"Datos MIBiG"/"mibig_gbk_4.0"
mibig_diamond = base/"Datos MIBiG"/"mibig_prot_seqs_4.0.dmnd"

json_anti = (base/"GCA_028473745.1_ASM2847374v1_genomic.json")

arch_reg = entrada_reg = resultados /"processed_regions.csv"

#Read output from 01_validate_regions

def leer_regionesmibig():
    if not arch_reg.exists():
        print("Error: The file does not exist")
        print(arch_reg)
        sys.exit(1)

    regiones = []

    with open(arch_reg, encoding="utf-8", newline="") as archivo:
        lector = csv.DictReader(archivo)

        if "Region" not in lector.fieldnames:
            print("Error: The CSV file does not contain the 'Region' column")
            sys.exit(1)

        if "Next Step" not in lector.fieldnames:
            print("Error: The CSV file does not contain the 'Next Step' column")  
            sys.exit(1)

        for fila in lector:
            sgte = fila["Next Step"].strip()

            if sgte == "MIBiG":
                regiones.append(fila["Region"].strip()) 
    return regiones

#Upload antiSMASH JSON
def cargar_json(ruta):
    if not ruta.exists():
        print("ERROR: The antiSMASH JSON does not exist:")
        print(ruta)
        sys.exit(1)

    with open(ruta, "r", encoding="utf-8") as archivo:
        return json.load(archivo)

#Get information about the region
def obtener_region_anti(data, region):
    partes = region.split(".")
    record_id = ".".join(partes[:-1])
    region_nombre = partes[-1]

    try:
        region_numero = int(region_nombre.replace("region", ""))
    except ValueError:
        return None, None, None

    for record in data.get("records", []):

        if record.get("id") != record_id:
            continue

        areas = record.get("areas", [])

        indice = region_numero - 1

        if indice < 0 or indice >= len(areas):
            return record, None, None

        area = areas[indice]

        modules = record.get("modules", {})

        cluster_compare = modules.get("antismash.modules.cluster_compare", {})

        db_results = cluster_compare.get("db_results", {})

        mibig = db_results.get("MIBiG", {})

        by_region = mibig.get("by_region",{})

        region_data = by_region.get(str(region_numero),{})

        return record, area, region_data

    return None, None, None

#Get core genes
def obtener_genecore(record, area):
    core_genes = []
    features = record.get("features", [])

    if not area:
        return core_genes

    protoclusters = area.get("protoclusters", {})

    if not protoclusters:
        return core_genes

    primer_protocluster = next(iter(protoclusters.values()))

    core_inicio = primer_protocluster.get("core_start")
    core_fin = primer_protocluster.get("core_end")

    if core_inicio is None or core_fin is None:
        return core_genes

    for feature in features:
        if not isinstance(feature, dict):
            continue

        if feature.get("type") != "CDS":
            continue

        qualifiers = feature.get("qualifiers", {})
        gene_kind = qualifiers.get("gene_kind", [])

        if "biosynthetic" not in gene_kind:
            continue

        location = feature.get("location", "")

        if not isinstance(location, str):
            continue

        match = re.search(r"\[(\d+):(\d+)\]", location)

        if not match:
            continue

        feature_inicio = int(match.group(1))
        feature_fin = int(match.group(2))

        if feature_inicio < core_inicio or feature_fin > core_fin:
            continue

        locus_tag = qualifiers.get("locus_tag", feature.get("id", ""))

        if isinstance(locus_tag, list):
            if locus_tag:
                locus_tag = locus_tag[0]
            else:
                locus_tag = ""

        if locus_tag:
            core_genes.append(str(locus_tag))

    return list(dict.fromkeys(core_genes))

#Region to Region
def obtener_r2r(region_data):
    resultados = []

    comparacion = region_data.get(
        "RegionToRegion_RiQ",
        {}
    )

    if not comparacion:
        return resultados

    scores = comparacion.get(
        "scores_by_region",
        {}
    )

    hits = comparacion.get(
        "hits",
        {}
    )

    for ref_id, score_bgc in scores.items():

        hits_bgc = hits.get(
            ref_id,
            {}
        )

        # The BGC does not contain any protein hits
        
        if not hits_bgc:

            resultados.append({
                "mibig_bgc": ref_id.split(":")[0],
                "score": score_bgc,
                "hits": []
            })

            continue

        # A BGC can have multiple protein hits

        hits_lista = []

        for proteina, hit in hits_bgc.items():

            hits_lista.append({
                "proteina_mibig": proteina,
                "query_cds_name": hit.get(
                    "query_cds_name",
                    ""
                ),
                "percent_identity": hit.get(
                    "percent_identity",
                    ""
                ),
                "percent_coverage": hit.get(
                    "percent_coverage",
                    ""
                ),
                "evalue": hit.get(
                    "evalue",
                    ""
                ),
                "blast_score": hit.get(
                    "blast_score",
                    ""
                )
            })

        resultados.append({
            "mibig_bgc": ref_id.split(":")[0],
            "score": score_bgc,
            "hits": hits_lista
        })

    #Sort by score from highest to lowest
    resultados.sort(
        key=lambda x: float(
            x["score"]
        ) if x["score"] not in [
            None,
            ""
        ] else 0,
        reverse=True
    )

    return resultados

#Extract the genes supported by each MIBiG BGC
def obtener_genes_supp(hits):
    genes = []
    for hit in hits:
        qseqid = hit.get("query_cds_name", "")
        if qseqid:
            genes.append(qseqid)
    return list(dict.fromkeys(genes))
     

#Read DIAMOND results
def leer_diamond(ruta):
    resultados = []

    columnas = [
        "qseqid",
        "sseqid",
        "pident",
        "length",
        "mismatch",
        "gapopen",
        "qstart",
        "qend",
        "sstart",
        "send",
        "evalue",
        "bitscore",
        "qcovhsp",
        "scovhsp"
    ]

    if not ruta.exists():
        return resultados

    with open(ruta,"r",encoding="utf-8") as archivo:

        lector = csv.reader(archivo,delimiter="\t")

        for fila in lector:

            if not fila:
                continue

            registro = dict(zip(columnas, fila))

            resultados.append(registro)

    return resultados

#Run DIAMOND
def ejecutar_diamond(region):
    fasta = (base / "results" / "MIBiG Homology" / f"{region}_proteinas.fasta")

    dir_salida_diamond = (base / "results" / "MIBiG Homology")

    dir_salida_diamond.mkdir( parents=True,exist_ok=True)

    archivo_tsv = (dir_salida_diamond/ f"{region}_homologia.tsv")

    if not fasta.exists():

        print(f"WARNING: There is no FASTA for"f"{region}:")
        print(f"  {fasta}")

        return []

    if not mibig_diamond.exists():

        print("WARNING: There is no MIBiG DIAMOND database:")
        print(f"  {mibig_diamond}")

        return []

    comando = [
        "diamond",
        "blastp",

        "--query",
        str(fasta),

        "--db",
        str(mibig_diamond),

        "--out",
        str(archivo_tsv),

        "--outfmt",
        "6",
        "qseqid",
        "sseqid",
        "pident",
        "length",
        "mismatch",
        "gapopen",
        "qstart",
        "qend",
        "sstart",
        "send",
        "evalue",
        "bitscore",
        "qcovhsp",
        "scovhsp",

        "--evalue",
        "1e-5",

        "--id",
        "50",

        "--query-cover",
        "70",

        "--max-target-seqs",
        "20"
    ]

    print(f"\nRunning DIAMOND for {region}...")

    resultado = subprocess.run(
        comando,
        capture_output=True,
        text=True
    )

    if resultado.returncode != 0:

        print("ERROR while running DIAMOND:")
        print(resultado.stderr)

        return []

    return leer_diamond(archivo_tsv)

#Group DIAMOND results by BGC MIBiG
def agrupar_diamond_bgc(hits):
    agrupados = {}

    for hit in hits:

        sseqid = hit.get("sseqid", "")

        if not sseqid:
            continue

        # MIBiG identifiers typically contain
        # the BGC accession number and the protein identifier.
        mibig_bgc = sseqid.split("|")[0]

        if mibig_bgc not in agrupados:
            agrupados[mibig_bgc] = []

        agrupados[mibig_bgc].append(hit)

    return agrupados

# Read information from the MIBiG JSON
def cargar_mibigjson(accession):
    accession_base = accession.split(".")[0]
    archivo = mibig_json / f"{accession_base}.json"

    if not archivo.exists():
        print(f"WARNING: There is no JSON MIBiG for {accession}")
        return None

    try:
        with open(archivo, "r", encoding="utf-8") as f:
            return json.load(f)

    except Exception as e:
        print(f"WARNING: Could not read {archivo}: {e}")

        return None

#Extract information from the compound
def extraer_infomibig(data):

    if not data:
        return {
            "compuesto": "",
            "smiles": "",
            "formula": "",
            "peso_molecular": "",
            "bioactividad": ""
        }

    compounds = data.get("compounds", [])

    nombres = []
    smiles = []
    formulas = []
    pesos = []
    bioactividades = []

    for compound in compounds:

        #Compound name

        nombre = compound.get("name", "")

        if nombre:
            nombres.append(str(nombre))

        #SMILES structure

        structure = compound.get("structure", "")

        if isinstance(structure, str) and structure.strip():
            smiles.append(structure.strip())

        # Molecular formula

        formula = compound.get("formula", "")

        if formula:
            formulas.append(str(formula))

        #Molecular weigth

        mass = compound.get("mass", "")

        if mass not in [None, ""]:
            pesos.append(str(mass))

        #Bioactivity

        actividades = compound.get("bioactivities", [])

        for actividad in actividades:

            if isinstance(actividad, dict):

                activity_data = actividad.get("name", "")

                if isinstance(activity_data, dict):
                    texto = (
                        activity_data.get("activity")
                        or activity_data.get("description")
                        or ""
                    )

                else:
                    texto = str(activity_data)

                if texto:
                    bioactividades.append(texto)

    return {
        "compuesto": "; ".join(dict.fromkeys(nombres)),
        "smiles": "; ".join(dict.fromkeys(smiles)),
        "formula": "; ".join(dict.fromkeys(formulas)),
        "peso_molecular": "; ".join(dict.fromkeys(pesos)),
        "bioactividad": "; ".join(dict.fromkeys(bioactividades))
    }

#Operational decision
def decidir_priorizacion(score_r2r, estructura_disponible, n_core_genes, n_core_genes_soportados):
    #Pipeline operational rule.

    # STRUCTURAL_CANDIDATE:
        #   - RegionToRegion evidence exists
        #   - A structure is available in MIBiG
        #   - Core genes are defined
        #   - At least 50% of the core genes are supported
        #
        # MANUAL_REVIEW_CANDIDATE:
        #   - partial evidence exists, but not all criteria are met
        #
        # UNASSIGNED:
        #   - no RegionToRegion evidence exists

    # 1. RegionToRegion evidence must exist

    # RegionToRegion evidence must exist
    if score_r2r in [None, ""]:
        return "NOT_ASSIGNED"

    #There must be a structure in place at MIBiG
    if not estructura_disponible:
        return "CANDIDATE_MANUAL_REVIEW"

    # There must be defined core genes
    if n_core_genes == 0:
        return "CANDIDATE_MANUAL_REVIEW"

    #Percentage of supported core genes
    proporcion_core = (n_core_genes_soportados / n_core_genes)

    #  Operational criterion: at least 50%
    if proporcion_core >= 0.50:
        return "STRUCTURAL_CANDIDATE"

    return "CANDIDATE_MANUAL_REVIEW"
#Process a region
def procesar_region(region, data_antismash):
    print("\n" + "=" * 70)
    print(f"REGION: {region}")
    print("=" * 70)

    record, area, region_data = (obtener_region_anti(data_antismash, region))

    if record is None:
        print("ERROR: Region not found in the antiSMASH JSON")
        return [], []

    if area is None:
        print("ERROR: The corresponding area was not found")
        return [], []

    # Product

    productos = area.get("products", [])

    print(f"Product(s): {', '.join(productos)}")

    #Core genes

    core_genes = obtener_genecore(record, area)

    print( f"Core genes: {len(core_genes)}")

    print(f"  {core_genes}")

    # RegionToRegion

    candidatos_r2r = obtener_r2r(region_data)

    print(f"RegionToRegion: {len(candidatos_r2r)} candidatos")


    resultados = []
    estructuras = []

    # BGC Candidates from RegionToRegion

    for candidato in candidatos_r2r:

        mibig_bgc = candidato["mibig_bgc"]

        score = candidato["score"]

        hits_r2r = candidato["hits"]

        genes_r2r = obtener_genes_supp(hits_r2r)

        # MIBiG information

        mibig_data = cargar_mibigjson(mibig_bgc)

        info = extraer_infomibig(mibig_data)

        estructura_disponible = bool(info["smiles"])

        # Genes supported by R2R
        genes_soportados_r2r = list(dict.fromkeys(genes_r2r))

        # Distinction between core genes and additional genes
        genes_core_soportados_r2r = [
            gen for gen in genes_soportados_r2r
            if gen in core_genes
        ]


        
        # Final consolidation
        genes_core_soportados = list(dict.fromkeys(
            genes_core_soportados_r2r
        ))

   
        # Decision
        
        decision = decidir_priorizacion(
            score_r2r=score,
            estructura_disponible=estructura_disponible,
            n_core_genes=len(core_genes),
            n_core_genes_soportados=len(genes_core_soportados)
        )

        print(
            f"DEBUG | {region} | {mibig_bgc} | "
            f"R2R={score} | "
            f"SMILES={'Yes' if estructura_disponible else 'No'} | "
            f"Core={len(core_genes)} | "
            f"Supported cores ={len(genes_core_soportados)} | "
            f"Supported core genes ={genes_core_soportados} | "
        )

        # Registration

        fila = {

            "Region":
                region,

            "Region_product":
                "; ".join(productos),

            "Core_genes":
                "; ".join(core_genes),

            "BGC_MIBiG":
                mibig_bgc,

            "MIBiG_compound":
                info["compuesto"],

            "SMILES":
                info["smiles"],

            "Formula":
                info["formula"],

            "Molecular_weight":
                info["peso_molecular"],

            "Bioactivity":
                info["bioactividad"],

          
            "N_core_genes":
                len(core_genes),

            "Score_RegionToRegion":
                score,

            "Structure_available":
                "Yes"
                if estructura_disponible
                else "No",

            "Decision":
                decision
        }

        resultados.append(
            fila
        )


    #Selection of the prioritized putative structure 

    candidatos_estructurales = [
        fila for fila in resultados
        if fila["Decision"] == "STRUCTURAL_CANDIDATE"
    ]

    if candidatos_estructurales:

        # Prioritize the candidate with the most evidence
        # antiSMASH RegionToRegion
        candidatos_estructurales.sort(
            key=lambda x: float(x["Score_RegionToRegion"])
            if x["Score_RegionToRegion"] not in [None, ""]
            else 0,
            reverse=True
        )

        seleccionado = candidatos_estructurales[0].copy()

        # The selected candidate becomes the
        # prioritized putative structure for the region
        seleccionado["Decision"] = "PRIORITIZED_PUTATIVE"

        estructuras.append(seleccionado)
    return resultados, estructuras

#Save Results
def guardar_csv(archivo, datos):

    if not datos:
        print(f"\nThere are no results to save in:")
        print(archivo)
        return

    archivo.parent.mkdir(parents=True, exist_ok=True)

    campos = list(datos[0].keys())

    with open(archivo, "w", newline="", encoding="utf-8-sig") as f:

        writer = csv.DictWriter(f, fieldnames=campos)

        writer.writeheader()

        writer.writerows(datos)

    print(f"\nFile generated:")
    print(archivo)

#MAIN
if __name__ == "__main__":

    print("=" * 70)
    print("       STRUCTURAL ALLOCATION - PHASE II")
    print("=" * 70)

    # Read regions from script 01

    regiones = leer_regionesmibig()

    print(f"\nRegions with an MIBiG route: {len(regiones)}")

    for region in regiones:
        print(f"  - {region}")

    if not regiones:

        print("\nNo regions were found with 'Next step = MIBiG'.")
        sys.exit(0)

    # Load antiSMASH


    data_antismash = cargar_json(json_anti)

    # Process

    todos_resultados = []
    todas_estructuras = []

    for region in regiones:

        resultados, estructuras = (procesar_region(region, data_antismash))
        todos_resultados.extend(resultados)
        todas_estructuras.extend(estructuras)

        
    # Save

    guardar_csv(sal_asignacion,todos_resultados)

    guardar_csv(sal_estructuras, todas_estructuras)

    # Summary

    print("\n" + "=" * 70)
    print("                 SUMMARY")
    print("=" * 70)

    print(f"Regions processed: {len(regiones)}")
    print(f"Candidates evaluated: {len(todos_resultados)}")
    print(f"Prioritized putative structures: {len(todas_estructuras)}")
    print("\nProcess completed")
