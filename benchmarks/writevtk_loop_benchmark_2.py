import time
import numpy as np

elnames = ['A', 'B'] * 5
phnames = ['P1', 'P2'] * 5

mf = [np.random.rand(100, 100) for _ in elnames]
mur = [np.random.rand(100, 100) for _ in elnames]
phf = [np.random.rand(100, 100) for _ in phnames]

np.set_printoptions(threshold=np.inf, linewidth=100000)

def old_way():
    h_str = 'POINT_DATA 10000\n'
    k = 0
    for i in elnames:
        h_str = h_str + 'SCALARS ' + 'mole-fraction(' + str(i) + ') Double 1'+'\n'+'LOOKUP_TABLE default' + '\n'  \
                        + np.array_str(mf[k][:]).replace('[', '').replace(']', '')+'\n' \
                        + 'SCALARS ' + 'chemical-potential(' + str(i) + ') Double 1'+'\n'+'LOOKUP_TABLE default'+'\n' \
                        + np.array_str(mur[k][:]).replace('[', '').replace(']', '')+'\n'
        k += 1
    k = 0
    for i in phnames:
        h_str = h_str + 'SCALARS ' + 'phase-fraction(' + str(i) + ') Double 1'+'\n'+'LOOKUP_TABLE default' + '\n' \
                        + np.array_str(phf[k][:]).replace('[', '').replace(']', '')+'\n'
        k += 1
    return h_str

def new_way_array_str():
    parts = ['POINT_DATA 10000\n']
    for k, i in enumerate(elnames):
        parts.append('SCALARS mole-fraction(' + str(i) + ') Double 1\nLOOKUP_TABLE default\n' \
                     + np.array_str(mf[k][:]).replace('[', '').replace(']', '') + '\n' \
                     + 'SCALARS chemical-potential(' + str(i) + ') Double 1\nLOOKUP_TABLE default\n' \
                     + np.array_str(mur[k][:]).replace('[', '').replace(']', '') + '\n')
    for k, i in enumerate(phnames):
        parts.append('SCALARS phase-fraction(' + str(i) + ') Double 1\nLOOKUP_TABLE default\n' \
                     + np.array_str(phf[k][:]).replace('[', '').replace(']', '') + '\n')
    return "".join(parts)

def new_way_ravel():
    parts = ['POINT_DATA 10000\n']
    for k, i in enumerate(elnames):
        parts.append('SCALARS mole-fraction(' + str(i) + ') Double 1\nLOOKUP_TABLE default\n' \
                     + ' '.join(map(str, mf[k].ravel())) + '\n' \
                     + 'SCALARS chemical-potential(' + str(i) + ') Double 1\nLOOKUP_TABLE default\n' \
                     + ' '.join(map(str, mur[k].ravel())) + '\n')
    for k, i in enumerate(phnames):
        parts.append('SCALARS phase-fraction(' + str(i) + ') Double 1\nLOOKUP_TABLE default\n' \
                     + ' '.join(map(str, phf[k].ravel())) + '\n')
    return "".join(parts)


start = time.time()
r1 = old_way()
end1 = time.time() - start

start = time.time()
r2 = new_way_array_str()
end2 = time.time() - start

start = time.time()
r3 = new_way_ravel()
end3 = time.time() - start

print(f"Old way: {end1:.4f}s")
print(f"New way array_str: {end2:.4f}s")
print(f"New way ravel: {end3:.4f}s")
