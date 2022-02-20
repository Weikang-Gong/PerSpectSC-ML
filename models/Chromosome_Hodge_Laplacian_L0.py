import numpy as np
import gudhi as gd
from scipy.sparse import dok_matrix
from scipy.io import loadmat

def n_faces(face_set, n):
    return filter(lambda face: len(face)==n+1, face_set)
    
def boundary_operator(face_set, i):
    source_simplices = list(n_faces(face_set, i))
    target_simplices = list(n_faces(face_set, i-1))

    if len(target_simplices)==0:
        S = dok_matrix((1, len(source_simplices)), dtype=np.float64)
        S[0, 0:len(source_simplices)] = 1
    else:
        source_simplices_dict = {source_simplices[j]: j for j in range(len(source_simplices))}
        target_simplices_dict = {target_simplices[i]: i for i in range(len(target_simplices))}

        S = dok_matrix((len(target_simplices), len(source_simplices)), dtype=np.float64)
        for source_simplex in source_simplices:
            for a in range(len(source_simplex)):
                target_simplex = source_simplex[:a]+source_simplex[(a+1):]
                i = target_simplices_dict[target_simplex]
                j = source_simplices_dict[source_simplex]
                S[i, j] = -1 if a % 2==1 else 1
    return S

DM=loadmat('H1_ESC_Chromosome_Distance_Matrix_K.mat')
Chr_Num=len(DM['H1_ESC_Chromosome_Distance_Matrix_K'][0])
str1='H1_TB_Chromosome_'
str3='_VR_L0_EV.txt'
for Chr in range(Chr_Num):
    if (Chr<=21):
        str2=str(Chr+1)
    elif (Chr==22):
        str2="X"
    elif (Chr==23):
        str2="Y"
    filename=str1+str2+str3;
    file=open(filename,'w')
    for i in range(101):
        filtration_value=i/100;
        rc = gd.RipsComplex(distance_matrix=DM['H1_ESC_Chromosome_Distance_Matrix_K'][0][Chr], max_edge_length=filtration_value)
        simplex_tree = rc.create_simplex_tree(max_dimension=1)
        val = simplex_tree.get_filtration()
        simplices = set()
        for v in val:
            simplices.add(tuple(v[0]))
        laplacian = np.matmul(boundary_operator(simplices, 1).toarray(), np.transpose(boundary_operator(simplices, 1).toarray()))
        eigval, eigvec = np.linalg.eigh(laplacian)
        eigvalue=np.sort(eigval)
        filtration_value='%.2f'%filtration_value
        file.write(str(filtration_value)+'\t')
        n=len(eigvalue)
        for j in range(n):
            EV='%.3f'%eigvalue[j]
            file.write(str(EV)+'\t')
        file.write('\n')
    file.close()