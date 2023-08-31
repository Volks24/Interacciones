import requests
import pandas as pd
import os


### Descargar PDB desde archivo 'pdb_descarga.csv' 
# necesita una columna denominada 'Estructura' para obtener los PDB 

if __name__ == '__main__':
    
    PDB_DF = pd.read_csv('pdb_descarga.csv')
    
    PDB_list = PDB_DF['Estructura'].tolist()
    
    if not os.path.exists('PDB'):
        os.makedirs('PDB')

    for PDB in PDB_list:
     
        # URL del archivo PDB que deseas descargar
        pdb_url = "https://www.rcsb.org/pdb/download/downloadFile.do?fileFormat=pdb&compression=NO&structureId={}".format(PDB)


        # Nombre del archivo de salida
        output_file = "PDB/{}.pdb".format(PDB)

        # Realizar la solicitud HTTP y descargar el archivo
        response = requests.get(pdb_url)

        if response.status_code == 200:
            with open(output_file, 'wb') as f:
                f.write(response.content)
        else:
            print("Error al descargar el pdb {}".format(PDB))