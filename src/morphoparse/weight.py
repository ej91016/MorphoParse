from pathlib import Path

def write_paup_weights(output_file, paup_lines):
    paup_file = f"{Path(output_file)}_paup_weights.txt"
    with open(paup_file, 'w') as f:
        f.write("BEGIN PAUP;\n")
        for line in paup_lines:
            weight, rest = line.split(":")
            positions = rest.strip().split()
            current_line = f"  WEIGHTS {weight}:"
            for pos in positions:
                if len(current_line) + len(pos) + 1 > 250:
                    f.write(current_line + ";\n")
                    current_line = f"  WEIGHTS {weight}: {pos}"
                else:
                    current_line += f" {pos}"
            f.write(current_line + ";\n")
        f.write("END;\n")
    print(f"SSA weight scheme written to {Path(paup_file).name}")

def write_tnt_weights(output_file, tnt_lines):
    tnt_file = f"{Path(output_file)}_tnt_weights.txt"
    with open(tnt_file, 'w') as f:
        f.write("\n".join(tnt_lines) + "\n")
    print(f"SSA weight scheme written to {Path(tnt_file).name}")