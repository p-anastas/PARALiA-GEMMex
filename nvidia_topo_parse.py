import pandas as pd
import numpy as np
import subprocess
import os

system_name = os.environ.get("PARALIA_GEMMEX_SYSTEM")
if not system_name:
    raise SystemExit('PARALIA_GEMMEX_SYSTEM  not defined, fill and run config_system.sh prior to this script!') 
# Define the path
system_dir = os.path.join('Deployment_files', system_name)
# Create the directory
os.makedirs(system_dir, exist_ok=True)  # exist_ok=True prevents an error if the directory already exists

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
chl_workers = len(parsed_smi)

def get_unique_values(df, column_name):
    # Get unique values from the specified column and convert to a list
    unique_values = df[column_name].unique().tolist()
    return unique_values

def translate_unit_list_to_binary(active_unit_id_list, max_len):
    case_id_out = 0

    for mask_offset in range(max_len):
        if mask_offset in active_unit_id_list:
            case_id_out += 2 ** mask_offset

    return case_id_out

dev_list = get_unique_values(parsed_smi,'dev_id')
#print(dev_list)
dev_case_id = translate_unit_list_to_binary(dev_list, chl_workers)
print(dev_case_id)

with open('%s/chl_worker_grid_%d.log' %(system_dir, dev_case_id), "w") as f:
    print('====================================================================================================', file=f)
    
    print('CHL_WORKERS = %d\n' % chl_workers, file=f)

    chl_dtypes = len(parsed_smi['gpu_model_flops'][0][0])
    print('WORKER_GOPS: %d' % chl_dtypes, file=f)

    for i in range(chl_dtypes):
        name = parsed_smi['gpu_model_flops'][0][0][i]
        print('%s :' % name, end='', file=f)
        
        for j in range(chl_workers):
            if name != parsed_smi['gpu_model_flops'][j][0][i]:
                raise SystemExit('Parsing bug, %s != %s' % (name, parsed_smi['gpu_model_flops'][j][0][i]))
            
            ops = parsed_smi['gpu_model_flops'][j][1][i]
            print(' %d' % ops, end='', file=f)
        print(file=f)

    print('\nWORKER_POWER:\nWATTS:', end='', file=f)
    for j in range(chl_workers):
        wats = parsed_smi['max_power_w'][j]
        print(' %d' % wats, end='', file=f)
    print('\n\n====================================================================================================', file=f)

subprocess.run(['cat', '%s/chl_worker_grid_%d.log' %(system_dir, dev_case_id)])

def parse_nvbandwidth_data(filename,breakline,lines_read):
    nvbandwidth_data = []
    
    with open(filename, 'r') as file:
        lines = file.readlines()
        
        # Find the starting line of the matrix data section
        start_line = 0
        for i, line in enumerate(lines):
            if breakline in line:
                start_line = i + 2  # The data starts two lines after this one
                break
        
        # Read and process each matrix row from the file
        for i in range(start_line, start_line + lines_read):
            row_data = lines[i].split()[1:]  # Ignore the first element (row index)
            parsed_row = []
            for value in row_data:
                # Replace "N/A" with a specific value, e.g., -1
                if value == "N/A":
                    parsed_row.append(-1)
                else:
                    parsed_row.append(round(float(value),2))  # Convert to int after parsing float
            nvbandwidth_data.append(parsed_row)

    return nvbandwidth_data

def balance_float_list_contents(list_bws):
    new_list = []
    for sublist in list_bws:
        # Calculate the average of the list
        avg = round(sum(sublist) / len(sublist), 2)
        # Replace each element with the average
        new_list.append([avg] * len(sublist))
    return new_list

def colonize_list2_in_list_of_lists(list_of_lists, list2):
    # Determine how many lists need to be filled with '-1'
    num_lists = len(list_of_lists)
    
    # Extend each sublist in list_of_lists
    for i in range(num_lists):
        if i < len(list2):
            # Append the element from list2 in order
            list_of_lists[i].append(list2[i])
        else:
            # Append '-1' for remaining lists
            list_of_lists[i].append(-1)

    return list_of_lists

def homogenize_all_values(list_of_lists, tolerance=0.05):
    # Flatten the list of lists into a single list to work with all values at once
    flat_list = [item for sublist in list_of_lists for item in sublist]
    groups = []  # List to hold groups of values that are within the tolerance
    # Group values based on the tolerance
    for value in flat_list:
        added_to_group = False
        for group in groups:
            # Check if the current value differs by less than tolerance from the group's average
            if abs(value - np.average(group)) / abs(np.average(group)) <= tolerance:
                group.append(value)
                added_to_group = True
                break
        if not added_to_group:
            groups.append([value])  # Create a new group if no match is found
    #print(groups)
    # Calculate the average for each group
    averages = {value: sum(group) / len(group) for group in groups for value in group}
    
    # Replace all values in list_of_lists with their group average
    homogenized_list = [[round(averages[value],2) for value in sublist] for sublist in list_of_lists]
    
    return homogenized_list

# Load data
nvbandwidth_data_d2d = homogenize_all_values(parse_nvbandwidth_data(system_dir + '/microbench_logs/' +
    'device_to_device_memcpy_write_ce.log', breakline='memcpy CE GPU(row) <- GPU(column) bandwidth (GB/s)',lines_read = chl_workers))
nvbandwidth_data_d2d_bid = homogenize_all_values(parse_nvbandwidth_data(system_dir + '/microbench_logs/' +
    'device_to_device_bidirectional_memcpy_write_ce.log', breakline='memcpy CE GPU(row) <-> GPU(column) Write1 bandwidth (GB/s)',lines_read = chl_workers))

nvbandwidth_data_d2h = balance_float_list_contents(parse_nvbandwidth_data(system_dir + '/microbench_logs/' +
    'device_to_host_memcpy_ce.log', breakline='memcpy CE CPU(row) <- GPU(column) bandwidth (GB/s)',lines_read = 1))
nvbandwidth_data_d2h_bid = balance_float_list_contents(parse_nvbandwidth_data(system_dir + '/microbench_logs/' +
    'device_to_host_bidirectional_memcpy_ce.log', breakline='memcpy CE CPU(row) <-> GPU(column) bandwidth (GB/s)',lines_read = 1))

nvbandwidth_data_d2h_inter = balance_float_list_contents(parse_nvbandwidth_data(system_dir + '/microbench_logs/' +
    'device_to_host_memcpy_ce_inter.log', breakline='memcpy CE CPU(row) <- GPU(column) bandwidth (GB/s)',lines_read = 1))
nvbandwidth_data_d2h_bid_inter = balance_float_list_contents(parse_nvbandwidth_data(system_dir + '/microbench_logs/' +
    'device_to_host_bidirectional_memcpy_ce_inter.log', breakline='memcpy CE CPU(row) <-> GPU(column) bandwidth (GB/s)',lines_read = 1))

nvbandwidth_data = nvbandwidth_data_d2d + nvbandwidth_data_d2h + nvbandwidth_data_d2h_inter
nvbandwidth_data_bid = nvbandwidth_data_d2d_bid + nvbandwidth_data_d2h_bid + nvbandwidth_data_d2h_bid_inter

nvbandwidth_data_h2d = balance_float_list_contents(parse_nvbandwidth_data(system_dir + '/microbench_logs/' +
    'host_to_device_memcpy_ce.log', breakline='memcpy CE CPU(row) -> GPU(column) bandwidth (GB/s)',lines_read = 1))
nvbandwidth_data_h2d_bid = balance_float_list_contents(parse_nvbandwidth_data(system_dir + '/microbench_logs/' +
    'host_to_device_bidirectional_memcpy_ce.log', breakline='memcpy CE CPU(row) <-> GPU(column) bandwidth (GB/s)',lines_read = 1))
colonize_list2_in_list_of_lists(nvbandwidth_data,nvbandwidth_data_h2d[0])
colonize_list2_in_list_of_lists(nvbandwidth_data_bid,nvbandwidth_data_h2d_bid[0])

nvbandwidth_data_h2d_inter = balance_float_list_contents(parse_nvbandwidth_data(system_dir + '/microbench_logs/' +
    'host_to_device_memcpy_ce_inter.log', breakline='memcpy CE CPU(row) -> GPU(column) bandwidth (GB/s)',lines_read = 1))
nvbandwidth_data_h2d_bid_inter = balance_float_list_contents(parse_nvbandwidth_data(system_dir + '/microbench_logs/' +
    'host_to_device_bidirectional_memcpy_ce_inter.log', breakline='memcpy CE CPU(row) <-> GPU(column) bandwidth (GB/s)',lines_read = 1))
colonize_list2_in_list_of_lists(nvbandwidth_data,nvbandwidth_data_h2d_inter[0])
colonize_list2_in_list_of_lists(nvbandwidth_data_bid,nvbandwidth_data_h2d_bid_inter[0])

# Print to verify
#for row in nvbandwidth_data:
#    print(row)

#for row in nvbandwidth_data_bid:
#    print(row)

chl_memlocs = chl_workers + 2
with open('%s/chl_bw_grid_%d_%d.log' %(system_dir, dev_case_id, dev_case_id), 'w') as f:
    # Print header lines
    print("====================================================================================================", file=f)
    print("CHL_MEMLOCS = %d" % chl_memlocs, file=f)
    print("", file=f)
    print("One-directional:", file=f)

    # Print each row for the "One-directional" section from nvbandwidth_data
    for row in nvbandwidth_data:
        print(" ".join(map(str, row)), file=f)

    print("", file=f)
    print("Bidirectional:", file=f)

    # Print each row for the "Bidirectional" section from nvbandwidth_data_bid
    for row in nvbandwidth_data_bid:
        print(" ".join(map(str, row)), file=f)

    print("", file=f)
    print("====================================================================================================", file=f)

subprocess.run(['cat', '%s/chl_bw_grid_%d_%d.log' %(system_dir, dev_case_id, dev_case_id)])

### Note: These are not used currently
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

#devices, matrix = load_nvidia_topo_matrix()

# Display the results
#print("Devices:", devices)
#print("Matrix:")
#for row in matrix:
#   print(row)

def translate_binary_to_unit_list(case_id, chl_workers):
    active_unit_id_list = []
    active_unit_num = 0

    for mask_offset in range(chl_workers):
        mask = 1 << mask_offset
        if case_id & mask:
            active_unit_id_list.append(mask_offset)
            active_unit_num += 1

    return active_unit_num, active_unit_id_list

def translate_unit_list_to_binary(active_unit_id_list, chl_workers):
    case_id_out = 0

    for mask_offset in range(chl_workers):
        if mask_offset in active_unit_id_list:
            case_id_out += 2 ** mask_offset

    return case_id_out

def binary_case_id_split(case_id, chl_workers):
    # Translate binary to list of active unit IDs
    active_unit_num, active_unit_id_list = translate_binary_to_unit_list(case_id, chl_workers)
    
    # Select every second active unit ID for the output list
    out_active_unit_id_list = active_unit_id_list[::2]  # Takes every second element
    
    # Check if the number of active units is even
    out_active_unit_num = len(out_active_unit_id_list)
    if out_active_unit_num * 2 != active_unit_num:
        raise ValueError(f"binary_case_id_split({case_id}): active_unit_num(={active_unit_num}) is not even")
    
    # Convert back to binary representation
    return translate_unit_list_to_binary(out_active_unit_id_list, chl_workers)

def is_subset(case_id, case_id_set, chl_workers):
    for mask_offset in range(chl_workers):
        mask = 1 << mask_offset
        # Check if case_id has a 1 where case_id_set has a 0
        if not (case_id_set & mask) and (case_id & mask):
            return False
    return True