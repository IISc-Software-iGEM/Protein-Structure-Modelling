def replace_coordinates(pdb_filename, new_coordinates_filename, output_filename):
    with open(pdb_filename, 'r') as pdb_file, open(new_coordinates_filename, 'r') as coord_file:
        pdb_lines = pdb_file.readlines()
        new_coords = coord_file.readlines()

    # Ensure we have enough new coordinates
    if len(new_coords) < len(pdb_lines):
        raise ValueError("Not enough new coordinates provided")

    with open(output_filename, 'w') as output_file:
        new_coord_index = 0
        for line in pdb_lines:
            if line.startswith('ATOM'):
                # Use new coordinates
                new_coord = new_coords[new_coord_index].strip().split()
                new_coord_index += 1
                # Format the new line
                new_line = f"{line[:30]}{float(new_coord[0]):8.3f}{float(new_coord[1]):8.3f}{float(new_coord[2]):8.3f}{line[54:]}"
                output_file.write(new_line + '\n')
            else:
                output_file.write(line)


# Example usage
pdb_filename = '../1qhk/input.pdb'
new_coordinates_filename = 'testfile.txt'
output_filename = 'test.pdb'

replace_coordinates(pdb_filename, new_coordinates_filename, output_filename)
