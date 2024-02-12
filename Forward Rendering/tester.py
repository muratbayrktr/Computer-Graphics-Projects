import os

def rasterizer(input_path, output_path):
    command = f"./code_template/rasterizer {input_path} -o {output_path}"
    os.system(command)

current_path = os.getcwd()

for root, dirs, files in os.walk(current_path):
    for file in files:
        if file.endswith(".xml"):
            input_path = os.path.join(root, file)

            # Create a new directory for the output
            output_dir = os.path.join(root, "rasterized_outputs")
            if not os.path.exists(output_dir):
                os.makedirs(output_dir)

            # Define the output file path
            output_file = os.path.splitext(file)[0] + "_output"  # Add your desired output file extension if needed
            output_path = os.path.join(output_dir, output_file)

            # Run the rasterizer
            rasterizer(input_path, output_path)
            print(f"Processed: {input_path}, Output: {output_path}")
     
            
            
