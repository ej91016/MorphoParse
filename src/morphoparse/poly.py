# Might expand if polymorphism get increased attention
class Poly:
    def __init__(self, poly_dict):
        self.data = poly_dict or {}

    def subset(self, selected_indices):
        index_map = {orig: new for new, orig in enumerate(selected_indices)}
        new_poly = {}

        for taxon, sites in self.data.items():
            subset_sites = {}
            for orig_index, states in sites.items():
                if orig_index in index_map:
                    subset_sites[index_map[orig_index]] = states
            if subset_sites:
                new_poly[taxon] = subset_sites

        return new_poly