import json
import JFIO

import numpy as np

def apply_pbc(x, y, z, Lx, Ly, Lz):
    """
    Aplica condiciones de contorno periódicas (PBC) a las coordenadas (x, y, z)
    dentro de una caja cúbica de tamaño L.

    Parámetros:
    x, y, z - Coordenadas del punto.
    L - Longitud del lado de la caja cúbica.

    Retorna:
    Coordenadas (x, y, z) ajustadas por las PBC.
    """
    def pbc(coord, L):
        return ((coord + L/2) % (L)) - L/2
    
    x_pbc = pbc(x, Lx)
    y_pbc = pbc(y, Ly)
    z_pbc = pbc(z, Lz)
    
    return x, y, z

with open("./parameters.json") as f:
    inputData = json.load(f)

tol = 1e-4

sigma  = inputData["sigma"]
epsilon = inputData["epsilon"]

vertices = inputData["vertex"]
faces = inputData["faces"]

N = inputData["N"]

Lx = inputData["Lx"]
Ly = inputData["Ly"]
Lz = inputData["Lz"]

T  = inputData["T"]
dt = inputData["dt"]

nSteps = inputData["nSteps"]

# Read the simulation data

sim = JFIO.read("./results/simulation.json")

# Read the simulation data
file_path = "./results/potential.dat"
data = np.loadtxt(file_path,skiprows=2)
E_sim = data[:,1]

def is_point_in_polygon(point, polygon):
    n = len(polygon)
    angle_sum = 0
    
    for i in range(n):
        v1 = polygon[i] - point
        v2 = polygon[(i + 1) % n] - point
        angle = np.arccos(np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2)))
        angle_sum += angle
    
    return np.isclose(angle_sum, 2 * np.pi)

def distance_point_to_edge(point, edge_start, edge_end):
    point = np.array(point)
    edge_start = np.array(edge_start)
    edge_end = np.array(edge_end)
    edge_vector = edge_end - edge_start
    edge_vector = edge_vector/np.linalg.norm(edge_vector)
    point_vector = point - edge_start
    gradient = 0
    
    t = np.dot(point_vector, edge_vector) / np.dot(edge_vector, edge_vector)
    
    if t < 0.0:
        nearest_point = edge_start

    elif t > 1.0:
        nearest_point = edge_end
        
    else:     
        nearest_point = edge_start + t * edge_vector
        
    dist = np.linalg.norm(point - nearest_point)
    gradient = (point - nearest_point)/dist
    
    return dist,gradient

def distance_point_to_face(point, face_vertices):
    gradient = 0
    point = np.array(point)
    vertices = np.array(face_vertices)
    
    normal = np.cross(vertices[1] - vertices[0], vertices[2] - vertices[0])
    normal = normal / np.linalg.norm(normal)
    
    point_to_face_vector = point-vertices[0]
    distance_to_plane = np.dot(point_to_face_vector, normal)
    projection = point - distance_to_plane * normal
    
    if is_point_in_polygon(projection, vertices):
        gradient = normal*np.sign(distance_to_plane)
        return np.abs(distance_to_plane),gradient
    
    min_distance = float('inf')
    for i in range(len(vertices)):
        edge_start = vertices[i]
        edge_end = vertices[(i + 1) % len(vertices)]
        aux_distance,aux_gradient = distance_point_to_edge(point, edge_start, edge_end)
        if min_distance > aux_distance:
            min_distance = aux_distance
            gradient = aux_gradient
    
    return min_distance, gradient

def minimum_distance_point_to_polyhedron(point, vertices, faces):
    min_distance = float('inf')
    gradient = 0
    
    # for vertex in vertices:
    #     min_distance = min(min_distance, distance_point_to_vertex(point, vertex))
    
    # for edge in edges:
    #     edge_start, edge_end = vertices[edge[0]], vertices[edge[1]]
    #     min_distance = min(min_distance, distance_point_to_edge(point, edge_start, edge_end))
    
    for face in faces:
        face_vertices = [vertices[i] for i in face]
        aux_distance,aux_gradient = distance_point_to_face(point, face_vertices)
        if min_distance > aux_distance:
            min_distance = aux_distance
            gradient = aux_gradient
    
    return min_distance,gradient

def gradientDistance(x,y,z):

    xp,yp,zp = apply_pbc(x,y,z,Lx,Ly,Lz)

    r_wall,_ = minimum_distance_point_to_polyhedron(np.array([xp,yp,zp]),vertices,faces)

    return r_wall*r_wall

def energy(x,y,z):

    gradDistance = gradientDistance(x,y,z);
    r2wall = 0
    r2wall = gradDistance*1

    if(r2wall > sigma*sigma*1.259921*1.259921):
        return 0.0;

    if(r2wall == 0.0):
        return np.nan

    invr2wall = sigma*sigma/r2wall;
    invr6wall = invr2wall*invr2wall*invr2wall;

    e = 4.0*epsilon*(invr6wall*invr6wall-invr6wall)+epsilon;

    return e



someError = False
state = sim["state"]["data"]
for particle in state:
    i = int(particle[0])

    x = float(particle[1][0])
    y = float(particle[1][1])
    z = float(particle[1][2])

    E_theo = energy(x,y,z)

    diff = np.abs((E_theo-E_sim[i])/(E_theo+E_sim[i]))
    if (diff > tol):
        someError = True
        print("Error for particle ",particle[0]," at  tion ",x,y,z)
        print("Theoretical energy: ",E_theo," Simulation energy: ",E_sim[i])
        print("Error: ",diff)
        print("")

if not someError:
    print("No errors found in the simulation")



