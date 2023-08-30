### prepare pdb ### 
### split pdb file for analisys ###

## Librerias ##
import argparse
import sys

## Funciones ##

def armar_pdb_out(PDB_Input ,pdb_file,others_pdb):

    PDB_OUT = open('r_{}.pdb'.format(pdb_file), 'w')
    
    for lines in PDB_Input:
        if 'ATOM' in (lines[0:4]):
            PDB_OUT.write(lines)
    PDB_OUT.write('TER\n')
    if len(others_pdb) > 0:
        for other in others_pdb:
            for lines in PDB_Input:
                if other in (lines[17:20]):
                    PDB_OUT.write(lines)
            PDB_OUT.write('TER\n')
    PDB_OUT.close()

def armar_ligando_out(PDB_Input,ligand_file_not_H):
    
    ligand_file_not_H = ligand_file_not_H.upper()
    
    ligand_out = open('{}.pdb'.format(ligand_file_not_H), 'w')
    
    for lines in PDB_Input:
        if ligand_file_not_H in (lines[17:20]):
            ligand_out.write(lines)
    ligand_out.close()


    
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='genera pdb de receptor y ligando desde pdb')

    # Agregar argumentos
    parser.add_argument('-p', '--pdb_file', type=str, help='PDB File')
    parser.add_argument('-o', '--others', type=str, help='Others')
    parser.add_argument('-l', '--ligand_pdb', type=str, help='Ligand name')

    # Parsear los argumentos de la línea de comandos
    args = parser.parse_args()

    if not (args.pdb_file and args.others and args.ligand_pdb):
        parser.print_help()
        sys.exit(1)

    # Acceder a las variables ingresadas
    pdb_file = args.pdb_file
    others_pdb = (args.others).split(',')
    ligand_file_not_H = args.ligand_pdb

    ## Open PDB ##

    PDB_Input =open(pdb_file).readlines()
    
    armar_pdb_out(PDB_Input,pdb_file.split('.')[0],others_pdb)

    armar_ligando_out(PDB_Input,ligand_file_not_H)

