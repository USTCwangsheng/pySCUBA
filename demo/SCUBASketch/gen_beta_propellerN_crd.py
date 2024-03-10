import os, sys
import numpy as np
  
  
def to_cartesian(polar_vector):
    length, angle = polar_vector[0], polar_vector[1]
    return np.array([length*np.cos(angle), length*np.sin(angle)]).T


def gen_alphabeta_sketch(motif_num=7, beta_dist=5, beta_to_beta_dist= 3, beta_sheet_num=4, beta_length=5, N_pad=3, N_pad_minus_z=4, out_f=None):
    pseudo_motif_num = motif_num+1
    motif_angle = 360.0/pseudo_motif_num
    side0_angle = motif_angle/2
    side1_angle = (180 - motif_angle)/2
    side2_length = np.sin(np.deg2rad(side1_angle))/np.sin(np.deg2rad(side0_angle)) * beta_dist/2
    structure_rad = np.sqrt(side2_length ** 2 + (beta_dist/2) ** 2)
    
    motif_angle = np.linspace(0, 2 * np.pi, pseudo_motif_num)
    beta_motif_len = np.ones(pseudo_motif_num) * structure_rad
    beta_xy_coords = to_cartesian([beta_motif_len, motif_angle])[:-1]
    
    for sheet_idx in range(beta_sheet_num-1):
        beta_motif_len = np.ones(pseudo_motif_num) * (structure_rad + beta_to_beta_dist * (sheet_idx+1))
        new_beta_xy_coords = to_cartesian([beta_motif_len, motif_angle])[:-1]
        beta_xy_coords = np.concatenate([beta_xy_coords, new_beta_xy_coords], -1)
    
    xycoords = beta_xy_coords.reshape(-1, 2)
    xyzcoords = np.concatenate([xycoords, np.zeros(xycoords.shape[0])[:, None]], -1)
        
    with open(out_f, 'w') as writer:
        for motif_idx, motif in enumerate(xyzcoords):
            if ((motif_idx % 2) == 0):
                ss = 'E'
                ss_length = beta_length
                if (motif_idx == 0):
                    ss_length = beta_length + N_pad
                    z_coords = -N_pad_minus_z
                else:
                    z_coords = 0.0
                writer.write(f'{ss} {ss_length} {ss_length} N;{round(motif[0], 2)} {round(motif[1], 2)} {z_coords}; 0 0 -1\n') #### could not be easily changed, coupled with clock-wise 
            else:
                ss = 'E'
                ss_length = beta_length
            
                writer.write(f'{ss} {ss_length} {ss_length} C;{round(motif[0], 2)} {round(motif[1], 2)} {round(motif[2], 2)}; 0 0 1\n')
    
    
if __name__ == '__main__':    
    gen_dir = './beta_propellerN'

    cur_dir = os.getcwd()
    for motif_num in np.arange(9, 12, 1): #motif num range:8,9,10
        for beta_length in np.arange(4, 6, 1): #beta_length range:4,5residues
            os.makedirs(f'{gen_dir}/{motif_num}_{beta_length}', exist_ok=True)
            out_f = f'{gen_dir}/{motif_num}_{beta_length}/{motif_num}_{beta_length}.txt'
            sketch_par_f = f'{gen_dir}/{motif_num}_{beta_length}/sketch.par'
            with open(sketch_par_f, 'w') as writer:
                writer.write(f'START SketchPar\n')
                writer.write(f'SketchFile = {motif_num}_{beta_length}.txt\n')
                writer.write(f'LinkLoop = 1\n')
                writer.write(f'RandomSeed = 123\n')
                writer.write(f'OutputFile = {motif_num}_{beta_length}_\n')
                writer.write(f'GenerationNumber = 3\n')
                writer.write(f'END SketchPar\n')

            gen_alphabeta_sketch(motif_num, beta_dist=5, beta_to_beta_dist=5, beta_length=beta_length, out_f=out_f)
            os.chdir(f'{gen_dir}/{motif_num}_{beta_length}')
            os.system(f'~/workspace/pySCUBA/pySCUBA/cpp_bin/SCUBASketch sketch.par') #generate pdb from crd. Change this path to your installation path.
            os.chdir(cur_dir)
                
