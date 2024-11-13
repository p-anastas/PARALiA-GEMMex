import subprocess
import pandas as pd

def nvidia_gpu_gflops(model_string):
    if model_string.startswith('A100'):
        return [('MM_FP64','MM_FP32','MM_FP16'),(19500,19500,19500)]
    else:
        raise SystemExit('Unknown GPU model flops for %s' %(model_string))

def nvidia_smi_parser():
    # Run the nvidia-smi command and capture its output
    result = subprocess.run(['nvidia-smi'], capture_output=True, text=True)
    
    # Check if the command ran successfully
    if result.returncode != 0:
        print("Error running nvidia-smi:", result.stderr)
        return []
    
    # Split output into lines
    lines = result.stdout.splitlines()
    
    # List to hold information for each GPU
    gpu_data = []
    next_line_candidate = 0
    
    # Define column names
    column_names = ["dev_id", "gpu_model", "gpu_model_flops", "max_power_w", "max_mem_gb"]

    # Initialize an empty DataFrame with the specified columns
    df = pd.DataFrame(columns=column_names)

    # Parse each line for GPU ID, Name, Memory, and Max Power Cap
    for line in lines:
        if line == '+-----------------------------------------+------------------------+----------------------+' or \
            line == '|=========================================+========================+======================|':
            next_line_candidate = 1
            continue
        if next_line_candidate == 1:
            if line.startswith('|'):
                #print(line)
                #Dev id
                gpu_data.append(int(line.split()[1]))
                #GPU model
                gpu_data.append(line.split()[3])
                gpu_data.append(nvidia_gpu_gflops(gpu_data[-1]))
                next_line_candidate = 2
                continue
            else:
                next_line_candidate = 0
        if next_line_candidate == 2:
                #print(line + '\n')
                #Max powa
                gpu_data.append(int(line.split()[6][:-1]))
                #Max mem
                gpu_data.append(int(int(line.split()[10][:-3])/1024))
                #print (gpu_data)
                new_row = pd.DataFrame([gpu_data], columns=column_names)
                df = pd.concat([df, new_row], ignore_index=True)
                gpu_data = []
                next_line_candidate = 0 
    return df

# Fetch and print GPU information
parsed_smi = nvidia_smi_parser()
print(parsed_smi)

print('====================================================================================================')
chl_workers = len(parsed_smi)
print('CHL_WORKERS = %d\n' %(chl_workers))

chl_dtypes = len(parsed_smi['gpu_model_flops'][0][0])
print('WORKER_GOPS = %d' %(chl_dtypes))

for i in range(0,chl_dtypes):
    name = parsed_smi['gpu_model_flops'][0][0][i]
    print ('%s :' %name, end='')
    for j in range(0,chl_workers):
        if name != parsed_smi['gpu_model_flops'][j][0][i]:
            raise SystemExit('Parsing bug, %s != %s' %(name, parsed_smi['gpu_model_flops'][j][0][i]))
        ops = parsed_smi['gpu_model_flops'][j][1][i]
        print (' %d' %ops, end='')
    print()
print('\nWORKER_POWER:\nWATTS:', end='')
for j in range(0,chl_workers):
    wats = parsed_smi['max_power_w'][j]
    print (' %d' %wats, end='')
print('\n\n====================================================================================================')


def load_nvidia_topo_matrix():
    # Run the nvidia-smi topo --matrix command and capture the output
    result = subprocess.run(['nvidia-smi', 'topo', '--matrix'], capture_output=True, text=True)
    
    # Check if the command ran successfully
    if result.returncode != 0:
        print("Error running nvidia-smi:", result.stderr)
        return None, None
    
    # Split the output into lines
    lines = result.stdout.splitlines()

    # Remove header information and extract GPU matrix part
    matrix = []
    devices = []
    
    for line in lines:
        # Skip empty lines
        if not line.strip():
            continue
        
        # Start reading the matrix after the header row
        if line.startswith("GPU"):
            headers = line.strip().split()
            
    return devices, matrix

devices, matrix = load_nvidia_topo_matrix()

# Display the results
print("Devices:", devices)
print("Matrix:")
for row in matrix:
    print(row)

print('====================================================================================================')
print('CHL_MEMLOCS = %d\n' %(len(devices) + 2))

print('\n====================================================================================================')
