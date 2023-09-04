### prepare pdb ### 
### split pdb file for analisys ###

## Librerias ##
import argparse
import sys
import os

## Funciones ##

def armar_pdb_out(PDB_Input ,pdb_file,others_pdb,cadena_interes):

    PDB_OUT = open('{}'.format(pdb_file), 'w')
    
    for lines in PDB_Input:
        if ('ATOM' in (lines[0:4])) and (cadena_interes == lines[21] ) :
            PDB_OUT.write(lines)
    #PDB_OUT.write('TER\n')
    if len(others_pdb) > 0:
        for other in others_pdb:
            for lines in PDB_Input:
                if 'ANISOU' in lines:
                    pass
                elif (other in (lines[17:20])) and (cadena_interes == lines[21] ):
                    PDB_OUT.write(lines)
    #        PDB_OUT.write('TER\n')
    PDB_OUT.close()

def armar_ligando_out(PDB_Input,ligand_file_not_H,cadena_interes):
    
    ligand_file_not_H = ligand_file_not_H.upper()
    
    ligand_out = open('{}.pdb'.format(ligand_file_not_H), 'w')
    
    for lines in PDB_Input:
        if 'ANISOU' in lines:
            pass
        elif (ligand_file_not_H in (lines[17:20])) and (cadena_interes == (lines[21:22])):
            ligand_out.write(lines)
    ligand_out.close()



if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='genera pdb de receptor y ligando desde pdb')

    # Agregar argumentos
    parser.add_argument('-p', '--pdb_file', type=str, help='PDB File')
    parser.add_argument('-o', '--others', type=str, help='Others')
    parser.add_argument('-l', '--ligand_pdb', type=str, help='Ligand name')
    parser.add_argument('-c', '--cadena', type=str, help='cadena')

    # Parsear los argumentos de la línea de comandos
    args = parser.parse_args()

    if not (args.pdb_file and args.others and args.ligand_pdb):
        parser.print_help()
        sys.exit(1)

    # Acceder a las variables ingresadas
    pdb_file = args.pdb_file
    others_pdb = (args.others).split(',')
    ligand_file_not_H = args.ligand_pdb
    cadena_interes = args.cadena

    ## Open PDB ##
    pdb_file_original = pdb_file.split('.')[0]+'_All.pdb'
    os.rename(pdb_file,pdb_file_original)

    PDB_Input =open(pdb_file_original).readlines()
    
    armar_pdb_out(PDB_Input,pdb_file,others_pdb,cadena_interes)

    armar_ligando_out(PDB_Input,ligand_file_not_H,cadena_interes)

