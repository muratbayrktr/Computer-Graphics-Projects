import os
import numpy as np

# Run the render_all.sh file
# Measure the time
import time
start_time = time.time()
# os.system("./render_all.sh")
end_time = time.time()
print(f"Total Render time: {end_time - start_time:.2f} seconds")

# Define the acceptable error rate
error_rate = 1e-5

# Get the list of ppm files in the current directory
ppm_files = [f for f in os.listdir("./") if f.endswith(".ppm")]
# Sort
ppm_files.sort()
max_filename_length = max(len(f) for f in ppm_files)

# Create a dictionary to store the scores for each ppm file
scores = {}

# Compare the ppm files in the current directory with the corresponding files in ../outputs_nobfc/
for ppm_file in ppm_files:
    expected_file = "../outputs_nobfc/" + ppm_file
    num_diff_pixels = 0
    total_pixels = 0
    with open(ppm_file, "r") as f1, open(expected_file, "r") as f2:
        lines1 = f1.readlines()
        lines2 = f2.readlines()
        # Check first 3 lines are the same
        # First line is P3
        # Second line is width height
        # Third line is 255

        for i in range(3):
            if lines1[i] != lines2[i]:
                print(f"{ppm_file}: first 3 lines do not match")
                break

        width, height = map(int, lines1[1].split())
        max_pixel_value = int(lines1[2])

        pixels1 = np.array([list(map(int, line.split())) for line in lines1[3:]])
        pixels2 = np.array([list(map(int, line.split())) for line in lines2[3:]])

        diff_pixels = np.abs(pixels1 - pixels2)
        num_diff_pixels = np.sum(diff_pixels > error_rate * max_pixel_value)
        total_pixels = width * height

        # Calculate the score for this ppm file
        score = 1 - (num_diff_pixels / total_pixels)
        scores[ppm_file] = round(score * 100,4)
        
        if scores[ppm_file] >= 90:
            print(f"{ppm_file.ljust(max_filename_length)}:\t\033[92m{scores[ppm_file]:.2f}%\033[0m")
        else:
            print(f"{ppm_file.ljust(max_filename_length)}:\t\033[91m{scores[ppm_file]:.2f}%\033[0m")

# Calculate the final score by averaging the scores in the dictionary
final_score = round(sum(scores.values())/ len(scores),4) if len(scores) > 0 else 0

# Print the final score
if final_score >= 90:
    if len(scores) > 0:
        print("\n+{:-^25}+{:-^25}+".format("", ""))
        print("|{:^25}|{:^25}|".format("File Name", "Score"))
        print("+{:-^25}+{:-^25}+".format("", ""))
        for ppm_file, score in scores.items():
            print("|{:<25}|{:>25.2f}%|".format(ppm_file, score))
        print("+{:-^25}+{:-^25}+".format("", ""))
        print("|{:^25}|{:^25.2f}%|".format("Final Score", final_score))
        print("+{:-^25}+{:-^25}+".format("", ""))
    else:
        print("\nNo ppm files found in the current directory.")
else:
    print(f"\nFinal score: \033[91m{final_score:.2f}%\033[0m")
