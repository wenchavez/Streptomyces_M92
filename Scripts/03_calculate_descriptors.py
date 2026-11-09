from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski, Crippen, rdMolDescriptors
from rdkit.Chem import rdFingerprintGenerator
from pathlib import Path
import csv


#PATHs

ruta_script = Path(__file__).resolve().parent
ruta = ruta_script.parent
ruta_resultados = ruta /"results"


#  Structures obtained using PRISM
entrada_prism = ruta_resultados / "processed_regions.csv"

# Putative structures obtained by MIBiG homology modeling
entrada_mibig = (ruta_resultados/ "MIBiG Homology"/ "putative structures.csv")

#  Final File
output = ruta_resultados / "molecular_descriptors.csv"


# Parameters for Morgan Fingerprints
morgan_r = 2
morgan_nbits = 2048

fp_gen = rdFingerprintGenerator.GetMorganGenerator(radius=morgan_r, fpSize=morgan_nbits, includeChirality=False, useBondTypes=True, includeRingMembership=True)

# Validation and Canonicalization of SMILES

def validar_y_canonicalizar_smiles(smiles):
    """
    Validates a SMILES string using RDKit and, if it is valid,
    generates its canonical representation
    """

    if not smiles:
        return None

    smiles = smiles.strip()

    # Validation
    mol = Chem.MolFromSmiles(smiles)

    if mol is None:
        return None

    # Canonicalization
    smiles_canonical = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)

    return {"mol": mol, "SMILES_canonical": smiles_canonical}

# Calculation of Descriptors


def calcular_descriptores(mol):

    return {
        "Formula": rdMolDescriptors.CalcMolFormula(mol),
        "Molecular_weight": round(Descriptors.MolWt(mol), 3),
        "LogP": round(Crippen.MolLogP(mol), 3),
        "TPSA": round(rdMolDescriptors.CalcTPSA(mol), 3),
        "HBD": Lipinski.NumHDonors(mol),
        "HBA": Lipinski.NumHAcceptors(mol),
        "Rotatable_bonds": Lipinski.NumRotatableBonds(mol),
        "Rings": Lipinski.RingCount(mol)
    }


# Calculation of Morgan Fingerprints
def calcular_morgan_fingerprint(mol):

    fp = fp_gen.GetFingerprint(mol)

    return fp.ToBitString()

# Read PRISM structures

def leer_estructuras_prism():

    estructuras = []

    if not entrada_prism.exists():
        print("\nWARNING:")
        print(f"There is no file named PRISM: {entrada_prism}")
        return estructuras

    with open(entrada_prism, encoding="utf-8", newline="") as f:

        reader = csv.DictReader(f)

        for fila in reader:

            if fila.get("Status", "").strip() != "Complete":
                continue

            if fila.get("Next Step", "").strip() != "RDKit":
                continue

            smiles = fila.get("SMILES", "").strip()

            if not smiles:
                continue

            estructuras.append({
                "Region": fila.get("Region", "").strip(),
                "Source": "PRISM",
                "Scaffolds": fila.get("Scaffolds", "").strip(),
                "BGC_MIBiG": "",
                "Compound": "",
                "SMILES": smiles
            })

    return estructuras

# Read MIBiG structures

def leer_estructuras_mibig():

    estructuras = []

    if not entrada_mibig.exists():
        print("\nWARNING:")
        print(f"The MIBiG file does not exist: {entrada_mibig}")
        return estructuras

    with open(entrada_mibig, encoding="utf-8-sig", newline="") as f:

        reader = csv.DictReader(f)

        for fila in reader:

            region = fila.get("Region", "").strip()
            bgc = fila.get("BGC_MIBiG", "").strip()

            compuestos = fila.get("MIBiG_compound", "").strip()

            smiles = fila.get("SMILES", "").strip()

            if not smiles:
                continue

            lista_compuestos = [x.strip() for x in compuestos.split(";") if x.strip()]

            lista_smiles = [x.strip() for x in smiles.split(";") if x.strip()]

            if len(lista_compuestos) == len(lista_smiles):

                for compuesto, smi in zip(
                    lista_compuestos,
                    lista_smiles
                ):

                    estructuras.append({
                        "Region": region,
                        "Source": "MIBiG Homology",
                        "Scaffolds": "",
                        "BGC_MIBiG": bgc,
                        "Compound": compuesto,
                        "SMILES": smi
                    })

            else:

                for smi in lista_smiles:

                    estructuras.append({
                        "Region": region,
                        "Source": "MIBiG Homology",
                        "Scaffolds": "",
                        "BGC_MIBiG": bgc,
                        "Compound": compuestos,
                        "SMILES": smi
                    })

    return estructuras


# Process structures

def procesar_estructuras():

    estructuras_prism = leer_estructuras_prism()
    estructuras_mibig = leer_estructuras_mibig()

    estructuras = (estructuras_prism + estructuras_mibig)

    resultados = []

    print("\n" + "=" * 70)
    print("VALIDATION, CANONICALIZATION, AND DESCRIPTORS")
    print("=" * 70)

    print(f"\nPRISM structures: {len(estructuras_prism)}")
    print(f"MIBiG Homology Structures: {len(estructuras_mibig)}")
    print(f"Total number of structures: {len(estructuras)}")

    # Process each structure

    for i, fila in enumerate(estructuras, start=1):

        region = fila["Region"]
        fuente = fila["Source"]
        scaffold = fila["Scaffolds"]
        bgc = fila["BGC_MIBiG"]
        compuesto = fila["Compound"]
        smiles = fila["SMILES"]

        print("\n" + "-" * 70)
        print(f"Structure {i}/{len(estructuras)}")
        print(f"Region: {region}")
        print(f"Source: {fuente}")

        if scaffold:
            print(f"Scaffold: {scaffold}")

        if bgc:
            print(f"BGC MIBiG: {bgc}")

        if compuesto:
            print(f"Compound: {compuesto}")

        print(f"SMILES original: {smiles}")

        # Validation + canonicalization
       
        estructura = validar_y_canonicalizar_smiles(smiles)

        if estructura is None:

            print("Invalid SMILES")
            continue

        mol = estructura["mol"]
        smiles_canonical = estructura[
            "SMILES_canonical"
        ]

        print("Valid SMILES")
        print(
            f"SMILES canonical: "
            f"{smiles_canonical}"
        )

        #Descriptors
        
        descriptores = calcular_descriptores(mol)

        print(f"Formula: " f"{descriptores['Formula']}")

        print(f"Molecular weight: " f"{descriptores['Molecular_weight']:.3f}")

        print(f"LogP: " f"{descriptores['LogP']:.3f}")

        print(
            f"TPSA: " f"{descriptores['TPSA']:.3f}")

        print(f"HBD: " f"{descriptores['HBD']}")

        print(f"HBA: " f"{descriptores['HBA']}")

        print(f"Rotatable bonds: " f"{descriptores['Rotatable_bonds']}")

        print(f"Rings: " f"{descriptores['Rings']}")


        # Morgan fingerprints
        morgan = calcular_morgan_fingerprint(mol)

        # Save result

        resultado = {
            "Region": region,
            "Source": fuente,
            "Scaffolds": scaffold,
            "BGC_MIBiG": bgc,
            "Compound": compuesto,
            "SMILES_original": smiles,
            "SMILES_canonical": smiles_canonical,
            **descriptores,
            "Morgan_fingerprint": morgan
        }

        resultados.append(resultado)

    return resultados

# Save results

def guardar_resultados(resultados):

    campos = [
        "Region",
        "Source",
        "Scaffolds",
        "BGC_MIBiG",
        "Compound",
        "SMILES_original",
        "SMILES_canonical",
        "Formula",
        "Molecular_weight",
        "LogP",
        "TPSA",
        "HBD",
        "HBA",
        "Rotatable_bonds",
        "Rings",
        "Morgan_fingerprint"
    ] 
    output.parent.mkdir(parents=True, exist_ok=True)
    with open(output,"w",newline="",encoding="utf-8-sig") as f:

        writer = csv.DictWriter(f, fieldnames=campos)

        writer.writeheader()
        writer.writerows(resultados)

    print("\n" + "=" * 70)
    print("RESULTS SAVED")
    print("=" * 70)
    print(output)


# MAIN

if __name__ == "__main__":

    resultados = procesar_estructuras()

    guardar_resultados(resultados)

    print("\nProcess completed")