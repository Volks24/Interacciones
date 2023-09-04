

### Librerias ###
import dataframe_image as dfi
from rdkit import Chem
import pandas as pd
import numpy as np
import math
import yaml
from Bio.PDB import *
from Bio.SVDSuperimposer import SVDSuperimposer
import csv
import re
import argparse
import os
import sys

### Funciones ###

def carga_variables():
    # Cargo Variables Generales #
    with open(r'Interacciones_variables.yml') as file:
        Interaciones = yaml.load(file, Loader=yaml.FullLoader)

    Distances_Hidrogen_Bonds =float(Interaciones['distancias']['Distances_Hidrogen_Bonds'])
    Distances_Aromatic = float(Interaciones['distancias']['Distances_Aromatic'])
    Distancia_Hidrofobica = float(Interaciones['distancias']['Distances_Hidrofobica'])
    
    Aceptores_Prot = Interaciones['acceptors']
    Dadores_Prot = Interaciones['donors']
    Aceptot_antecedent = Interaciones['acceptors_antecedent']
    Special_case = Interaciones['special']
    return(Distances_Hidrogen_Bonds,Distances_Aromatic,Distancia_Hidrofobica,Aceptores_Prot,Dadores_Prot,Aceptot_antecedent,Special_case)

def Ligand_dataframe(ligando_estructura,ligand_name):
    ligand = pd.DataFrame(columns=['Ligand','ID' , 'Atom', 'X' , 'Y' , 'Z'])
    j = 0
    Coordenadas_interes = []
    for lines in ligando_estructura:
        if ('HETATM' in lines) or ('ATOM' in lines):
            ligand.loc[j , 'Ligand'] = ligand_name
            ligand.loc[j , 'Atom'] = str(lines[12:16])
            ligand.loc[j , 'ID'] = str(lines[7:11])
            ligand.loc[j , 'X'] = float(lines[31:38])
            ligand.loc[j , 'Y'] = float(lines[40:46])
            ligand.loc[j , 'Z'] = float(lines[47:55])
            Coordenadas_interes.append([float(lines[31:38]),float(lines[40:46]),float(lines[47:55])])
            j = j + 1
        else:
            pass
    center = calcular_centro_de_masa(Coordenadas_interes)
     
    return(ligand,center)

def calcular_centro_de_masa(puntos):
    sum_x = sum(punto[0] for punto in puntos)
    sum_y = sum(punto[1] for punto in puntos)
    sum_z = sum(punto[2] for punto in puntos)
    n = len(puntos)
    cm_x = sum_x / n
    cm_y = sum_y / n
    cm_z = sum_z / n
    return cm_x, cm_y, cm_z


def create_folders(carpetas):
    
    # Crear las carpetas una por una
    for carpeta in carpetas:
        if not os.path.exists(carpeta):  # Verificar si la carpeta no existe
            os.makedirs(carpeta)  # Crear la carpeta



def busqueda_antecesor(DF_Ligand , Pos ):
    ref_point = np.array(list(DF_Ligand.iloc[Pos,[3,4,5]])) # Punto de referencia
    sub_set = DF_Ligand.drop(DF_Ligand.index[Pos]) # lo saco para evitar el  0
    matriz_values = sub_set[['X', 'Y', 'Z']].values # matriz de datos
    # Calcular las distancias euclidianas entre el punto de referencia y cada punto en la matriz
    distancias = [math.sqrt((p[0]-ref_point[0])**2 + (p[1]-ref_point[1])**2 + (p[2]-ref_point[2])**2) for p in matriz_values]
    # Encontrar el índice del punto más cercano
    indice_punto_cercano = distancias.index(min(distancias))
    # Obtener el punto más cercano
    punto_cercano = matriz_values[indice_punto_cercano]
    Antecesor = (DF_Ligand.query('X == @punto_cercano[0]' and 'Y == @punto_cercano[1]' and 'Z == @punto_cercano[2]'))
    return(Antecesor)

def find_hidrogen_bond(DF_Ligand):
            
    ## Busco Candidatos a P.H. mediante 0 ##

    Posicion_Dador = []
    Posicion_Aceptor = []
    Antecesores_Aceptores = {}
    for Pos in range(0,DF_Ligand.shape[0]):  ## Busco oxigenos en ligando
        if 'O' in (DF_Ligand.iloc[Pos,2]):
            Cercano = busqueda_antecesor(DF_Ligand , Pos )
            Valor = (Cercano.iloc[0,2])
            if 'H' in Valor:
                Posicion_Dador.append(DF_Ligand.iloc[Pos,2].strip())
            Posicion_Aceptor.append(DF_Ligand.iloc[Pos,2].strip())
            Antecesores_Aceptores[DF_Ligand.iloc[Pos,2].strip()] = Valor.strip()
        if 'N' in (DF_Ligand.iloc[Pos,2]):
            Cercano = busqueda_antecesor(DF_Ligand , Pos )
            Valor = (Cercano.iloc[0,2])
            if 'H' in Valor:
                Posicion_Dador.append(DF_Ligand.iloc[Pos,2].strip())
            Posicion_Aceptor.append(DF_Ligand.iloc[Pos,2].strip())
            Antecesores_Aceptores[DF_Ligand.iloc[Pos,2].strip()] = Valor.strip()
        if 'F' in (DF_Ligand.iloc[Pos,2]):
            Cercano = busqueda_antecesor(DF_Ligand , Pos )
            Valor = (Cercano.iloc[0,2])
            if 'H' in Valor:
                Posicion_Dador.append(DF_Ligand.iloc[Pos,2].strip())
            Posicion_Aceptor.append(DF_Ligand.iloc[Pos,2].strip())
            Antecesores_Aceptores[DF_Ligand.iloc[Pos,2].strip()] = Valor.strip()
    #print('Aceptores {} \n Ante {} \n Dadores {}'.format(Posicion_Aceptor,Antecesores_Aceptores,Posicion_Dador))
    return(Posicion_Aceptor,Antecesores_Aceptores,Posicion_Dador)


def find_aromatic(mol,DF_Ligand):
    
    # Get the aromatic rings
    sssr = Chem.GetSSSR(mol)

    aromatic_rings = [ring for ring in sssr if len(ring) == 6]  # Filter out rings with fewer than 6 atoms
    # Obtener los átomos participantes de los anillos
    participating_rings = []

    for ring in aromatic_rings:
        participating_rings.append(list(ring))
    
 
    aromaticos = []

    for j in range(0,len(participating_rings)):
        r = []
        for k in range(0,len(participating_rings[j])):
            aro = (DF_Ligand.iloc[participating_rings[j][k],2])
            if any(char.isdigit() for char in aro):
                r.append(aro.strip())
            else:
                r.append(str(aro)+str((DF_Ligand.iloc[pos,1]).strip()))
        aromaticos.append(r)
    
    return(aromaticos)

def Normalize_DF_Ligand(DF_Ligand):
    Normalizar = 0
    Atomos = DF_Ligand['Atom'].tolist()
    for j in range(0,len(Atomos)):
        if 'H' in (Atomos[j]):
            Atomos[j] = Atomos[j].strip()+str(j+1)
        elif not re.search(r'\d', Atomos[j].strip()):
            Normalizar = 1
        else:
            Atomos[j] = Atomos[j].strip()
    if Normalizar == 1:
        DF_Ligand['ID_Atom'] = DF_Ligand.apply(lambda row: str(row['Atom']).strip() + str(row['ID']).strip(), axis=1)
        DF_Ligand['Atom'] = DF_Ligand['ID_Atom']
    else:
        DF_Ligand['Atom'] = Atomos
    return(DF_Ligand)


def output_file(Ligand_Name,Posicion_Aceptor,Posicion_Dador,Antecesores_Aceptores,Anillo):
    anillos_list = []
    anillos_par_list =  []
    anillos_par =  []
    archivo = open('{}/{}.yml'.format('YML',Ligand_Name) , 'w')
    archivo.write("### Datos Configuracion Interacciones\n\n")
    archivo.write("ligand:\n")
    archivo.write("  Atomos_Dadores : {}\n".format(Posicion_Dador))
    archivo.write("  Atomos_Aceptores: {}\n".format(Posicion_Aceptor))
    archivo.write("  Antecesores_Aceptores: {}\n".format(Antecesores_Aceptores))
    for j in range(0,len(Anillo)):
        anillos_list.append(Anillo[j])
        anillos_par_list.append([Anillo[j][0],Anillo[j][1]])
        anillos_par.append(Anillo[j][0])
        anillos_par.append(Anillo[j][1])
    archivo.write("  Atomos_Aromaticos: {} # 1 Pos x Anillo\n".format(anillos_list))
    archivo.write("  Anillo_Antecesor_Pares: {}\n".format(anillos_par_list))
    archivo.write("  Anillo_Antecesor: {}\n".format(anillos_par))
    archivo.write("  Especial: {}\n")
    archivo.close()


def active_site_residues(structure, Ligando_Centro,cadena, centroid_distance):
    model = structure[0][cadena]

    active_site = pd.DataFrame(columns=['Serial', 'Pos', 'Residue', 'Atom', 'X' , 'Y' ,'Z', 'CM X' , 'CM Y' , 'CM Z' ])

    Residuos_Interes = []

    for residue in model.get_residues():
        Residuo_Center = list(center_of_mass(residue))
        if (math.dist(Ligando_Centro, Residuo_Center)) < centroid_distance:
            Residuos_Interes.append([residue.get_resname(), residue.get_id()[1]])
            for atom in residue:
                Res_name = residue.get_resname()
                Res_id = residue.get_id()[1]
                atom_name = atom.get_name()
                Coor = list(atom.get_coord())
                Serial = atom.get_serial_number()
                #atoms.append([Serial, Res_id, Res_name, atom_name, Coor, Residuo_Center])
                active_site.loc[len(active_site.index)] = [Serial ,Res_id, Res_name, atom_name ,round(float(Coor[0]),3),round(float(Coor[1]),3),round(float(Coor[2]),3),round(Residuo_Center[0],3),round(Residuo_Center[1],3),round(Residuo_Center[2],3)]
    
    return active_site

def center_of_mass(entity, geometric=False):
    """
    Returns gravitic [default] or geometric center of mass of an Entity.
    Geometric assumes all masses are equal (geometric=True)
    """
    
    # Structure, Model, Chain, Residue
    if isinstance(entity, Entity.Entity):
        atom_list = entity.get_atoms()
    # List of Atoms
    elif hasattr(entity, '__iter__') and [x for x in entity if x.level == 'A']:
        atom_list = entity
    else: # Some other weirdo object
        raise ValueError("Center of Mass can only be calculated from the following objects:\n"
                            "Structure, Model, Chain, Residue, list of Atoms.")
    
    masses = []
    positions = [ [], [], [] ] # [ [X1, X2, ..] , [Y1, Y2, ...] , [Z1, Z2, ...] ]
    
    for atom in atom_list:
        masses.append(atom.mass)
        
        for i, coord in enumerate(atom.coord.tolist()):
            positions[i].append(coord)

    # If there is a single atom with undefined mass complain loudly.
    if 'ukn' in set(masses) and not geometric:
        raise ValueError("Some Atoms don't have an element assigned.\n"
                         "Try adding them manually or calculate the geometrical center of mass instead.")
    
    if geometric:
        return [sum(coord_list)/len(masses) for coord_list in positions]
    else:       
        w_pos = [ [], [], [] ]
        for atom_index, atom_mass in enumerate(masses):
            w_pos[0].append(positions[0][atom_index]*atom_mass)
            w_pos[1].append(positions[1][atom_index]*atom_mass)
            w_pos[2].append(positions[2][atom_index]*atom_mass)

        return [sum(coord_list)/sum(masses) for coord_list in w_pos]
    
def table_to_jpg(DF , out_put_name):
    
    if DF.shape[0] < 100:
                  
        df_styled = DF.style.background_gradient()
        dfi.export(df_styled, 'Graficos/{}.jpg'.format(out_put_name))
    else:
        # Define el máximo de filas por grupo
        max_rows_per_group = 100
        # Calcula el número total de grupos necesarios
        num_groups = (len(DF) + max_rows_per_group - 1) // max_rows_per_group

        # Divide el DataFrame en grupos más pequeños
        
        for i in range(num_groups):
            start_idx = i * max_rows_per_group
            end_idx = (i + 1) * max_rows_per_group
            df_group = DF.iloc[start_idx:end_idx, :]
            df_group = df_group.style.background_gradient()
            dfi.export(df_group, 'Graficos/{}_{}.jpg'.format(out_put_name,i))


    return()

def Coordenadas_interes_receptor(Aceptores_Prot,Dadores_Prot,DF_Active_Site):
    No_Enlista = ['HEM'] # Atomos no cubiertos
    ### Obtengo las coordenadas de los atomos de interes en el receptor
    receptor_points = pd.DataFrame(columns=['Type','Pos','Residue', 'Atom', 'X' , 'Y' , 'Z'])
    for pos in range(0,DF_Active_Site.shape[0]):
        Atomo = (DF_Active_Site.iloc[pos,2])
        if Atomo in No_Enlista:
            pass
        else:
            Res = (DF_Active_Site.iloc[pos,3])
            listado = (Aceptores_Prot[Atomo])
        if Res in listado:
            receptor_points.loc[len(receptor_points.index)] = 'Aceptor',DF_Active_Site.iloc[pos,1],DF_Active_Site.iloc[pos,2],DF_Active_Site.iloc[pos,3],DF_Active_Site.iloc[pos,4],DF_Active_Site.iloc[pos,5],DF_Active_Site.iloc[pos,6]
    for pos in range(0,DF_Active_Site.shape[0]):
        Atomo = (DF_Active_Site.iloc[pos,2])
        if Atomo in No_Enlista:
            pass
        else:
            Res = (DF_Active_Site.iloc[pos,3])
            listado = (Dadores_Prot[Atomo])
        if Res in listado:
            receptor_points.loc[len(receptor_points.index)] = 'Dador',DF_Active_Site.iloc[pos,1],DF_Active_Site.iloc[pos,2],DF_Active_Site.iloc[pos,3],DF_Active_Site.iloc[pos,4],DF_Active_Site.iloc[pos,5],DF_Active_Site.iloc[pos,6]
    aa_aro = ['TYR' , 'PHE' , 'TRP']
    for pos in range(0,DF_Active_Site.shape[0]):
        Atomo = (DF_Active_Site.iloc[pos,2])
        ID = (DF_Active_Site.iloc[pos,1])
        if Atomo in aa_aro:
            Sub_Set = DF_Active_Site.query('Pos == @ID')
            x,y,z = get_aromatic_coord(Atomo,Sub_Set)
            if ID not in receptor_points['Pos'].values:
                receptor_points.loc[len(receptor_points.index)] = 'Aromatico',DF_Active_Site.iloc[pos,1],DF_Active_Site.iloc[pos,2],'center',x,y,z
            elif 'Aromatico' not in (receptor_points.query('Pos == @ID')['Type'].tolist()) :# Solo posicion 
                receptor_points.loc[len(receptor_points.index)] = 'Aromatico',DF_Active_Site.iloc[pos,1],DF_Active_Site.iloc[pos,2],'center',x,y,z
        
            
    return(receptor_points)

def get_aromatic_coord(Res,AA):
    Aromatic_Ring = []
    if (Res == 'TYR') or (Res == 'PHE'):
            Coordenada = (AA.loc[AA['Atom'] == "CG", ['X','Y','Z']]).values.tolist()
            Aromatic_Ring.append(Coordenada[0])
            Coordenada = (AA.loc[AA['Atom'] == "CD1", ['X','Y','Z']]).values.tolist()
            Aromatic_Ring.append(Coordenada[0])
            Coordenada = (AA.loc[AA['Atom'] == "CD2", ['X','Y','Z']]).values.tolist()
            Aromatic_Ring.append(Coordenada[0])
            Coordenada = (AA.loc[AA['Atom'] == "CE1", ['X','Y','Z']]).values.tolist()
            Aromatic_Ring.append(Coordenada[0])
            Coordenada = (AA.loc[AA['Atom'] == "CE2", ['X','Y','Z']]).values.tolist()
            Aromatic_Ring.append(Coordenada[0])
            Coordenada = (AA.loc[AA['Atom'] == "CZ", ['X','Y','Z']]).values.tolist()
            Aromatic_Ring.append(Coordenada[0])
    elif ((Res) == 'TRP') :
            Coordenada = (AA.loc[AA['Atom'] == "CE3", ['X','Y','Z']]).values.tolist()
            Aromatic_Ring.append(Coordenada[0])
            Coordenada = (AA.loc[AA['Atom'] == "CD2", ['X','Y','Z']]).values.tolist()
            Aromatic_Ring.append(Coordenada[0])
            Coordenada = (AA.loc[AA['Atom'] == "CZ3", ['X','Y','Z']]).values.tolist()
            Aromatic_Ring.append(Coordenada[0])
            Coordenada = (AA.loc[AA['Atom'] == "CE2", ['X','Y','Z']]).values.tolist()
            Aromatic_Ring.append(Coordenada[0])
            Coordenada = (AA.loc[AA['Atom'] == "CH2", ['X','Y','Z']]).values.tolist()
            Aromatic_Ring.append(Coordenada[0])
            Coordenada = (AA.loc[AA['Atom'] == "CZ2", ['X','Y','Z']]).values.tolist()
            Aromatic_Ring.append(Coordenada[0])    
    center = center_aromatic_ring(Aromatic_Ring)
    return(center[0],center[1],center[2])
    

    
def center_aromatic_ring(Aromatic_Ring):
	x,y,z = [],[],[]

	for  j in range(0,len(Aromatic_Ring)):
		x.append(float(Aromatic_Ring[j][0]))
		y.append(float(Aromatic_Ring[j][1]))
		z.append(float(Aromatic_Ring[j][2]))
		
	CD1 = (x[1],y[1],z[1])
	CE1 = (x[3],y[3],z[3])

	vector_1 = (np.add(CD1, CE1))
	
	CD2 = (x[2],y[2],z[2])
	CE2 = (x[4],y[4],z[4])
		
	vector_2 = (np.add(CD2, CE2))
	
	center = (np.add(vector_2/2, vector_1/2))/2
	return(np.round(center,3))

def Distancia_Interaccion(Receptor , Ligando ,DF_Interacciones,DF_Ligand, Case , Distancia_PH,Antecesores_Aceptores):
    ref_point = np.array(list(Ligando.iloc[0,[3,4,5]])) # Punto de referencia
    matriz_values = Receptor[['X', 'Y', 'Z']].values # matriz de datos
    # Calcular las distancias euclidianas entre el punto de referencia y cada punto en la matriz
    distancias = [math.sqrt((p[0]-ref_point[0])**2 + (p[1]-ref_point[1])**2 + (p[2]-ref_point[2])**2) for p in matriz_values]
    for d in range(0,len(distancias)):
        if distancias[d] < Distancia_PH:
            if Case == 'Dador':
                Donor = list(Ligando.iloc[0,[3,4,5]]) # Ligando
                Aceptor = list(Receptor.iloc[d,[4,5,6]]) # Receptor
                Aceptor_Antecedent = list(Receptor.iloc[d-1,[4,5,6]]) # Receptor mas cercano p/angulo REVISAR QUE PUNTO LE DA
                angle = angle_three_points(Donor,Aceptor,Aceptor_Antecedent)
                if (angle > 100) and ((angle < 200)): ### Angulo optimo P.H.
                    Inter = 'Si'
                else:
                    Inter = 'No'
            if Case == 'Aceptor':
                Donor = list(Receptor.iloc[d,[4,5,6]]) # Receptor
                Aceptor = list(Ligando.iloc[0,[3,4,5]]) # Ligando
                atomo_ante = Antecesores_Aceptores[Ligando.iloc[0,2]]
                Aceptor_Antecedent = DF_Ligand.query('Atom == @atomo_ante')[['X' , 'Y' , 'Z']].values
                angle = angle_three_points(Donor,Aceptor,Aceptor_Antecedent[0])
                if (angle > 160) and ((angle < 200)): 
                    Inter = 'Si'
                else:
                    Inter = 'No' 
            DF_Interacciones.loc[len(DF_Interacciones.index)] = Receptor.iloc[d,1] , Receptor.iloc[d,2] ,Receptor.iloc[d,3] , distancias[d] , Ligando.iloc[0,2] , Case ,angle , Inter
    return(DF_Interacciones)        


def angle_three_points(Donor,Aceptor,Aceptor_Antecedent):
    
    Donor_coord = np.array(Donor)

    Aceptor_coord = np.array(Aceptor)

    Aceptor_Antecedent_coord = np.array(Aceptor_Antecedent)

    ba = Donor_coord-Aceptor_coord # normalization of vectors
    bc = Aceptor_Antecedent_coord-Aceptor_coord # normalization of vectors

    cosine_angle = np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc))
    angle = np.arccos(cosine_angle)
    
    return (np.degrees(angle))  # calculated angle in radians to degree 


def Interaccion_Aromatica(Receptor_Aromaticos,Aromaticos,DF_Interacciones,DF_Ligand,DF_Active_Site):
    for ring in Aromaticos:
        Anillo = (np.array(DF_Ligand[DF_Ligand['Atom'].isin(ring)][['X','Y','Z']]))
        Ligando_anillo_name = cadena = '-'.join(ring)
        Ligand_Center = center_aromatic_ring(Anillo)
        for k in range(0,Receptor_Aromaticos.shape[0]):
            Res = Receptor_Aromaticos.iloc[k,2]
            Pos = int(Receptor_Aromaticos.iloc[k,1]) # Pos ok       
            Sub_Set = DF_Active_Site[DF_Active_Site['Pos'] == Pos]
            if (Res == 'TYR') or (Res == 'PHE'): 
                Puntos_Interes = ['CG' , 'CD1' , 'CD2']
                Anillo_Name = 'CG-CD1-CD2'
                anillo_recept = np.array(Sub_Set[Sub_Set['Atom'].isin(Puntos_Interes)][['X','Y','Z']]).astype(float)
            elif (Res) == 'TRP':
                Puntos_Interes = ['CZ3' , 'CE3' , 'CH2']
                anillo_recept = np.array(Sub_Set[Sub_Set['Atom'].isin(Puntos_Interes)][['X','Y','Z']]).astype(float)
                Anillo_Name = 'CZ3-CE3-CH2'
            Aromatic_center_receptor = (np.array(Receptor_Aromaticos.iloc[k,[4,5,6]]))
            distance_centers = math.dist(Ligand_Center, Aromatic_center_receptor)
            if (distance_centers > 3) and (distance_centers < 4):
                angulo_planos = aromatic_angle(Anillo,anillo_recept)
                if (angulo_planos > 0) and (angulo_planos < 30):
                    Interaccion = 'Si'
                elif (angulo_planos > 85) and (angulo_planos < 95):
                    Interaccion = 'Si'
                else:
                    Interaccion = 'No'
                DF_Interacciones.loc[len(DF_Interacciones.index)]= Receptor_Aromaticos.iloc[k,1] , Receptor_Aromaticos.iloc[k,2] , Anillo_Name , distance_centers , Ligando_anillo_name,'Aromatico' , angulo_planos, Interaccion
    return(DF_Interacciones)
    
def aromatic_angle(Anillo,anillo_recept):
    anillo_ligand = np.array(Anillo[0:3]).astype(float)
    # Encontrar los átomos comunes más cercanos
    atomos_comunes = [anillo_ligand[0], anillo_recept[1]]
    # Calcular los vectores normales a los planos aromáticos
    vector_normal1 = np.cross(anillo_ligand[1] - atomos_comunes[0], anillo_ligand[2] - atomos_comunes[0])
    vector_normal2 = np.cross(anillo_recept[2] - atomos_comunes[1], anillo_recept[0] - atomos_comunes[1])
    # Calcular el ángulo entre los vectores normales
    producto_punto = np.dot(vector_normal1, vector_normal2)
    norma_vector1 = np.linalg.norm(vector_normal1)
    norma_vector2 = np.linalg.norm(vector_normal2)
    # Angulos 
    angulo_rad = np.arccos(producto_punto / (norma_vector1 * norma_vector2))
    angulo_deg = np.degrees(angulo_rad)
    return(angulo_deg)

def Interaccion_Especiales(DF_Interacciones,DF_Ligand,DF_Active_Site,special):
    
    Atom_especial = list(special.keys())
    for at in Atom_especial:
        Caso = (special[at])
        Recep = (DF_Active_Site.query('Atom == @Caso'))
        Point_rec = np.array(Recep[['X' , 'Y' , 'Z']])
        if len(Point_rec[0]) > 0:
            for llaves in Caso:
                resultado = DF_Ligand['Atom'].str.contains('N')
                Sub_Set = (DF_Ligand[resultado])
                for k in range(0,Sub_Set.shape[0]):
                    Point_lig = np.array(Sub_Set.iloc[k,[3,4,5]])
                    if (math.dist(Point_rec[0], Point_lig)) < float((special[at][1])):
                        DF_Interacciones.loc[len(DF_Interacciones.index)]= Recep.iloc[0,1] , Recep.iloc[0,2] , Recep.iloc[0,3] , math.dist(Point_rec[0], Point_lig) , Sub_Set.iloc[k,2] ,'Especial' , '-', 'Si'
    
    #DF_Interacciones = DF_Interacciones.drop((DF_Interacciones['Res'] == Recep.iloc[0,3]) and (DF_Interacciones['Interaccion'] == 'No'))
    DF_Interacciones = DF_Interacciones[~((DF_Interacciones['Res'] == 'HEM') & (DF_Interacciones['Interaccion'] == 'No'))]
    return(DF_Interacciones)




######### Funciones ###########

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Cargar pdb de receptor y')

    # Agregar argumentos
    parser.add_argument('-r', '--receptor_pdb', type=str, help='Receptor PDB')
    parser.add_argument('-l', '--ligand_pdb', type=str, help='Ligand PDB')

    # Parsear los argumentos de la línea de comandos
    args = parser.parse_args()

    if not (args.receptor_pdb and args.ligand_pdb):
        parser.print_help()
        sys.exit(1)

    # Acceder a las variables ingresadas
    receptor_pdb = args.receptor_pdb
    ligand_file_not_H = args.ligand_pdb

    # Protonate Ligand only polar #
    
    os.system(('obabel -ipdb {}.pdb -opdb -O {}_H.pdb -protonate polar').format(ligand_file_not_H.split('.')[0],ligand_file_not_H.split('.')[0]))
    
    ligand_file = ligand_file_not_H.split('.')[0]+'_H.pdb'

    Name_File = str(receptor_pdb.split('.')[0])+'_'+str(ligand_file.split('.')[0])

    ### Cargo datos externos ###

    Distances_Hidrogen_Bonds,Distances_Aromatic,Distancia_Hidrofobica,Aceptores_Prot,Dadores_Prot,Aceptot_antecedent,Special_case = carga_variables()
    Caso = ['acceptors','donors']
    create_folders(['YML' , 'Graficos' , 'Temp'])
    ## Crear Dataframe del Ligando ##

    ## Load the PDB file of ligand ##

    ligando_estructura = open('{}'.format(ligand_file),'r').readlines()
    ligand_case = ligand_file.split('.')[0]

    # Genero el DF #

    DF_Ligand , Centro_Masa = Ligand_dataframe(ligando_estructura,ligand_case)


    DF_Ligand = Normalize_DF_Ligand(DF_Ligand)


    # Calculo Hot points del ligando #

    Posicion_Aceptor,Antecesores_Aceptores,Posicion_Dador = find_hidrogen_bond(DF_Ligand)

    # Aromatic # 
    # Load the PDB file

    mol = Chem.MolFromPDBFile('{}'.format(ligand_file))
    
    Aromaticos = find_aromatic(mol,DF_Ligand) 
    print('{}\n{}\n{}\n{}'.format(Posicion_Aceptor,Antecesores_Aceptores,Posicion_Dador,Aromaticos))

    output_file(ligand_case,Posicion_Aceptor,Posicion_Dador,Antecesores_Aceptores,Aromaticos)

    ### Prepare Receptor ###

    ## Load the PDB file of receptor ##

    cadena_receptor = 'A'
    Distancia_Centro_Activo = 12


    pdb_parser = PDBParser()
    structure = pdb_parser.get_structure('pdb', receptor_pdb)

    ## ----------- Armo DataFrame Sitio Activo ----------- ##
    ## Obtengo dataframe con todo lo que este a X A° del centro de masa del ligando

    DF_Active_Site = active_site_residues(structure, Centro_Masa,cadena_receptor, Distancia_Centro_Activo)
    DF_Active_Site.to_csv('Temp/{}_Active_Site.csv'.format(str(receptor_pdb.split('.')[0])))

    ## Guardo como imagen ##
    table_to_jpg(DF_Active_Site.copy(),str(receptor_pdb.split('.')[0])+'_active_site')

    #### Interacciones #####
    
    receptor_points = Coordenadas_interes_receptor(Aceptores_Prot,Dadores_Prot,DF_Active_Site)

    ## ----------- Armo DataFrame Interacciones ----------- ##

    DF_Interacciones = pd.DataFrame( columns=['Pos R','Res','Atom' ,'Dist','Lig','P.H.', 'Angle' , 'Interaccion'])

    ### Calculos P.H: ###

    ## Ligando Dador - Receptor Aceptor ##

    Receptor_Caso ='Aceptor'
    Caso_Lig = 'Dador'

    Sub_Set = (receptor_points.query('Type == @Receptor_Caso'))

    for pos in Posicion_Dador:
        Ligand_res = (DF_Ligand.query('Atom == @pos'))
        DF_Interacciones = Distancia_Interaccion(DF_Active_Site , Ligand_res ,DF_Interacciones,DF_Ligand,Caso_Lig,Distances_Hidrogen_Bonds,Antecesores_Aceptores)


    ## Ligando Aceptor - Receptor Dador ##

    Receptor_Caso ='Dador'
    Caso_Lig = 'Aceptor'

    Sub_Set = (receptor_points.query('Type == @Receptor_Caso'))


    for pos in Posicion_Aceptor:
        Ligand_res = (DF_Ligand.query('Atom == @pos'))
        DF_Interacciones = Distancia_Interaccion(DF_Active_Site , Ligand_res ,DF_Interacciones,DF_Ligand,Caso_Lig,Distances_Hidrogen_Bonds,Antecesores_Aceptores)
    
    

    ## Aromaticas ##
    Receptor_Caso ='Aromatico'
    Receptor_Aromaticos = (receptor_points.query('Type == @Receptor_Caso'))
    DF_Interacciones  = Interaccion_Aromatica(Receptor_Aromaticos,Aromaticos,DF_Interacciones,DF_Ligand,DF_Active_Site)


    ## Especiales ##

    DF_Interacciones  = Interaccion_Especiales(DF_Interacciones,DF_Ligand,DF_Active_Site,Special_case)

    print(DF_Interacciones)

    ## Armo DF Final ##
    DF_Interacciones = DF_Interacciones.round(3)
    
    DF_Interacciones.to_csv('Temp/{}_{}.csv'.format(Name_File , 'Interaccion_All'))  

    Filtrado = (DF_Interacciones[DF_Interacciones['Interaccion'] == 'Si'])
    
    Filtrado.to_csv('Temp/{}_{}.csv'.format(Name_File , 'Interaccion'))  

    table_to_jpg(Filtrado.copy(),Name_File+'_interaccion')
    
    ##### archivo de salida resumen y chequear angulos 
    #### mejorar la parte de aromaticos query
