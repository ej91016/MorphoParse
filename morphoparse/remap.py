def remap_sparse_states(records, output_log_file, remove_missing=True, remove_uninformative=True, reorder_states=True):
    from Bio.SeqRecord import SeqRecord
    from Bio.Seq import Seq

    nchar = len(records[0].seq)
    char_buffers = [[] for _ in records]
    remapped_sites = []
    removed_sites = []

    for pos in range(nchar):
        col_chars = [r.seq[pos] for r in records]
        clean_chars = sorted({c for c in col_chars if c not in {'-', '?' }})
        invalid_chars = {c for c in clean_chars if not c.isalnum()}
        if invalid_chars:
            raise ValueError(f"Invalid characters found: {invalid_chars} at position {pos}")

        if (remove_missing and not clean_chars) or (remove_uninformative and len(clean_chars) == 1):
            removed_sites.append(pos + 1)
            continue

        state_map = {c: str(i) for i, c in enumerate(clean_chars)} if reorder_states else {c: c for c in clean_chars}
        if reorder_states and clean_chars != list("0123456789"[:len(clean_chars)]):
            remapped_sites.append((pos + 1, state_map.copy()))

        for i, r in enumerate(records):
            orig_char = col_chars[i]
            new_char = state_map.get(orig_char, orig_char)
            char_buffers[i].append(new_char)

    modified_records = [
        SeqRecord(Seq("".join(buf)), id=records[i].id)
        for i, buf in enumerate(char_buffers)
    ]

    with open(output_log_file, 'w') as f:
        for site, mapping in remapped_sites:
            mapping_str = ', '.join(f"{k}->{v}" for k, v in mapping.items())
            f.write(f"Site {site}: {mapping_str}\n")
        if removed_sites:
            f.write(f"\nRemoved {len(removed_sites)} site(s) with no data or single state: {', '.join(map(str, removed_sites))}\n")

    if remapped_sites or removed_sites:
        print(f"Remapping log written to: {output_log_file}")
    else:
        print("No remapping or column removal was necessary.")

    return modified_records