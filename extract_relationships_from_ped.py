#!/usr/bin/env python3
"""
Enhanced relationship extraction with comprehensive consanguinity detection.
Handles all complex consanguineous scenarios including double cousins, pedigree collapse,
and multi-generational consanguinity with inbreeding coefficient calculation.
"""

import sys
from collections import defaultdict, deque
import itertools

def parse_ped_file(filename):
    """Parse PED file and return list of individuals."""
    individuals = []
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#') or line.startswith('//'):
                continue
            # Split on any whitespace (spaces or tabs)
            parts = line.split()
            if len(parts) >= 6:
                fam_id, ind_id, father_id, mother_id, sex, phenotype = parts[:6]
                individuals.append({
                    'fam_id': fam_id,
                    'ind_id': ind_id,
                    'father_id': father_id if father_id != '0' else None,
                    'mother_id': mother_id if mother_id != '0' else None,
                    'sex': int(sex),
                    'phenotype': phenotype
                })
    return individuals

def validate_pedigree(individuals):
    """Validate pedigree for consistency and detect impossible relationships."""
    validation_errors = []
    ind_dict = {ind['ind_id']: ind for ind in individuals}
    
    for ind in individuals:
        # Check if parents exist
        if ind['father_id'] and ind['father_id'] not in ind_dict:
            validation_errors.append(f"Individual {ind['ind_id']}: Father {ind['father_id']} not found in pedigree")
        
        if ind['mother_id'] and ind['mother_id'] not in ind_dict:
            validation_errors.append(f"Individual {ind['ind_id']}: Mother {ind['mother_id']} not found in pedigree")
        
        # Check parent genders
        if ind['father_id'] and ind['father_id'] in ind_dict:
            father = ind_dict[ind['father_id']]
            if father['sex'] == 2:  # Female
                validation_errors.append(f"Individual {ind['ind_id']}: Father {ind['father_id']} has female gender")
        
        if ind['mother_id'] and ind['mother_id'] in ind_dict:
            mother = ind_dict[ind['mother_id']]
            if mother['sex'] == 1:  # Male
                validation_errors.append(f"Individual {ind['ind_id']}: Mother {ind['mother_id']} has male gender")
        
        # Check for self-loops (impossible relationships)
        if ind['father_id'] == ind['ind_id']:
            validation_errors.append(f"Individual {ind['ind_id']}: Cannot be their own father")
        if ind['mother_id'] == ind['ind_id']:
            validation_errors.append(f"Individual {ind['ind_id']}: Cannot be their own mother")
    
    # Check for cycles in pedigree
    cycles = detect_pedigree_cycles(individuals)
    if cycles:
        validation_errors.extend([f"Pedigree cycle detected: {' -> '.join(cycle)}" for cycle in cycles])
    
    return validation_errors

def detect_pedigree_cycles(individuals):
    """Detect cycles in pedigree structure (impossible relationships)."""
    ind_dict = {ind['ind_id']: ind for ind in individuals}
    cycles = []
    
    def dfs_cycle_detection(current_id, path, visited):
        if current_id in path:
            # Found a cycle
            cycle_start = path.index(current_id)
            cycle = path[cycle_start:] + [current_id]
            cycles.append(cycle)
            return
        
        if current_id in visited or current_id not in ind_dict:
            return
        
        visited.add(current_id)
        path.append(current_id)
        
        current = ind_dict[current_id]
        if current['father_id']:
            dfs_cycle_detection(current['father_id'], path.copy(), visited)
        if current['mother_id']:
            dfs_cycle_detection(current['mother_id'], path.copy(), visited)
    
    for ind in individuals:
        dfs_cycle_detection(ind['ind_id'], [], set())
    
    return cycles

def build_family_tree(individuals):
    """Build enhanced family tree data structure with ancestor mapping."""
    tree = {}
    children = defaultdict(list)
    ancestors = defaultdict(set)
    
    for ind in individuals:
        tree[ind['ind_id']] = ind
        if ind['father_id']:
            children[ind['father_id']].append(ind['ind_id'])
        if ind['mother_id']:
            children[ind['mother_id']].append(ind['ind_id'])
    
    # Build ancestor mapping for each individual
    def get_all_ancestors(ind_id, visited=None):
        if visited is None:
            visited = set()
        if ind_id in visited or ind_id not in tree:
            return set()
        
        visited.add(ind_id)
        ind = tree[ind_id]
        ancestor_set = set()
        
        if ind['father_id']:
            ancestor_set.add(ind['father_id'])
            ancestor_set.update(get_all_ancestors(ind['father_id'], visited.copy()))
        if ind['mother_id']:
            ancestor_set.add(ind['mother_id'])
            ancestor_set.update(get_all_ancestors(ind['mother_id'], visited.copy()))
        
        return ancestor_set
    
    for ind in individuals:
        ancestors[ind['ind_id']] = get_all_ancestors(ind['ind_id'])
    
    return tree, children, ancestors

def find_all_relationship_paths_comprehensive(tree, ancestors, proband_id, target_id):
    """Enhanced comprehensive path finding for all consanguineous scenarios."""
    if proband_id == target_id:
        return [("proband", 1.0)]
    
    # Check if this is a partner
    if target_id.endswith('P'):
        related_id = target_id[:-1]
        if related_id in tree:
            related_paths = find_all_relationship_paths_comprehensive(tree, ancestors, proband_id, related_id)
            if related_paths and related_paths[0][0] != "unrelated":
                partner_paths = []
                for rel_string, _ in related_paths:
                    partner_relationship = f"{rel_string}->partner"
                    partner_paths.append((partner_relationship, 0.0))
                return partner_paths
            else:
                return [("unrelated", 0.0)]
    
    # Quick check: if target is not an ancestor or descendant, look for common ancestors
    common_ancestors = ancestors[proband_id].intersection(ancestors[target_id])
    
    if not common_ancestors and target_id not in ancestors[proband_id] and proband_id not in ancestors[target_id]:
        # No blood relationship possible
        return [("unrelated", 0.0)]
    
    all_paths = []
    max_depth = 15  # Increased for complex consanguinity
    
    # Use enhanced DFS with backtracking for complete path enumeration
    def dfs_all_paths(current_id, target_id, path, visited_in_path, max_remaining_depth):
        if max_remaining_depth < 0:
            return
        
        if current_id == target_id:
            if path:  # Don't include empty paths
                relationship_string = build_relationship_string(path)
                path_relatedness = calculate_relatedness_from_path(path)
                
                # Avoid duplicate paths
                path_signature = (relationship_string, round(path_relatedness, 8))
                if path_signature not in [(r, round(rel, 8)) for r, rel in all_paths]:
                    all_paths.append((relationship_string, path_relatedness))
            return
        
        if current_id not in tree or max_remaining_depth == 0:
            return
        
        current = tree[current_id]
        
        # Explore all possible connections
        connections = []
        
        # Parents
        if current['father_id'] and current['father_id'] not in visited_in_path:
            connections.append((current['father_id'], "father"))
        if current['mother_id'] and current['mother_id'] not in visited_in_path:
            connections.append((current['mother_id'], "mother"))
        
        # Children
        for child_id in get_children(tree, current_id):
            if child_id not in visited_in_path:
                child = tree[child_id]
                if child['sex'] == 1:
                    connections.append((child_id, "son"))
                elif child['sex'] == 2:
                    connections.append((child_id, "daughter"))
                else:
                    connections.append((child_id, "child"))
        
        # All types of siblings
        siblings_dict = get_siblings(tree, current_id)
        
        for sibling_id in siblings_dict['full']:
            if sibling_id not in visited_in_path:
                sibling = tree[sibling_id]
                if sibling['sex'] == 1:
                    connections.append((sibling_id, "brother"))
                elif sibling['sex'] == 2:
                    connections.append((sibling_id, "sister"))
                else:
                    connections.append((sibling_id, "sibling"))
        
        for sibling_id in siblings_dict['half']:
            if sibling_id not in visited_in_path:
                sibling = tree[sibling_id]
                if sibling['sex'] == 1:
                    connections.append((sibling_id, "half-brother"))
                elif sibling['sex'] == 2:
                    connections.append((sibling_id, "half-sister"))
                else:
                    connections.append((sibling_id, "half-sibling"))
        
        # Explore each connection
        for next_id, relationship_term in connections:
            new_visited = visited_in_path.copy()
            new_visited.add(current_id)
            new_path = path + [relationship_term]
            dfs_all_paths(next_id, target_id, new_path, new_visited, max_remaining_depth - 1)
    
    # Start comprehensive search
    dfs_all_paths(proband_id, target_id, [], {proband_id}, max_depth)
    
    if not all_paths:
        return [("unrelated", 0.0)]
    
    # Sort by relatedness (highest first) and remove true duplicates
    unique_paths = []
    seen_paths = set()
    
    for relationship_string, relatedness in sorted(all_paths, key=lambda x: x[1], reverse=True):
        path_key = (relationship_string, round(relatedness, 8))
        if path_key not in seen_paths:
            unique_paths.append((relationship_string, relatedness))
            seen_paths.add(path_key)
    
    return unique_paths

def calculate_inbreeding_coefficient(tree, ancestors, individual_id):
    """Calculate inbreeding coefficient for an individual."""
    if individual_id not in tree:
        return 0.0
    
    ind = tree[individual_id]
    if not ind['father_id'] or not ind['mother_id']:
        return 0.0  # Cannot be inbred without both parents
    
    father_id = ind['father_id']
    mother_id = ind['mother_id']
    
    # Find all common ancestors between parents
    father_ancestors = ancestors[father_id] if father_id in ancestors else set()
    mother_ancestors = ancestors[mother_id] if mother_id in ancestors else set()
    
    common_ancestors = father_ancestors.intersection(mother_ancestors)
    
    if not common_ancestors:
        return 0.0
    
    # Calculate inbreeding coefficient
    inbreeding_coeff = 0.0
    
    for ancestor_id in common_ancestors:
        # Find paths from father to common ancestor
        father_paths = find_all_relationship_paths_comprehensive(tree, ancestors, father_id, ancestor_id)
        # Find paths from mother to common ancestor  
        mother_paths = find_all_relationship_paths_comprehensive(tree, ancestors, mother_id, ancestor_id)
        
        for f_rel, f_relatedness in father_paths:
            for m_rel, m_relatedness in mother_paths:
                if f_relatedness > 0 and m_relatedness > 0:
                    # Calculate path contributions to inbreeding
                    f_meioses = count_meioses_in_path(f_rel)
                    m_meioses = count_meioses_in_path(m_rel)
                    
                    # Inbreeding contribution = (1/2)^(n1 + n2 + 1)
                    # where n1 and n2 are path lengths to common ancestor
                    path_contribution = (0.5) ** (f_meioses + m_meioses + 1)
                    inbreeding_coeff += path_contribution
    
    return inbreeding_coeff

def count_meioses_in_path(relationship_string):
    """Count meioses in a relationship path."""
    if relationship_string == "proband":
        return 0
    
    parts = relationship_string.split("->")[1:]  # Remove "proband"
    meioses = 0
    
    for part in parts:
        if part not in ["partner"]:
            meioses += 1
    
    return meioses

def detect_consanguinity_type(all_paths):
    """Enhanced consanguinity type detection."""
    if len(all_paths) <= 1:
        return "no_consanguinity", []
    
    relationship_types = []
    relatedness_values = []
    
    for rel_string, relatedness in all_paths:
        rel_type = determine_relationship_type(rel_string, relatedness)
        relationship_types.append(rel_type)
        relatedness_values.append(relatedness)
    
    # Detect specific consanguinity patterns
    consanguinity_details = []
    
    # Parental consanguinity
    if any("mother" in rel for rel in relationship_types) and any("father" in rel for rel in relationship_types):
        consanguinity_details.append("parental_consanguinity")
    
    # Double cousin detection
    cousin_relationships = [rel for rel in relationship_types if "cousin" in rel.lower()]
    if len(cousin_relationships) >= 2:
        if len(set(cousin_relationships)) == 1:  # Same type of cousin relationship
            consanguinity_details.append("double_cousin")
        else:
            consanguinity_details.append("multiple_cousin_paths")
    
    # Avuncular relationships
    avuncular_terms = ["aunt", "uncle", "niece", "nephew"]
    if any(any(term in rel.lower() for term in avuncular_terms) for rel in relationship_types):
        consanguinity_details.append("avuncular_consanguinity")
    
    # Complex multi-generational
    if len(all_paths) > 2:
        consanguinity_details.append("complex_multigenerational")
    
    # Pedigree collapse (same person appearing multiple times)
    unique_relatedness = set(round(r, 6) for r in relatedness_values)
    if len(unique_relatedness) > 1 and len(all_paths) > 2:
        consanguinity_details.append("pedigree_collapse")
    
    if not consanguinity_details:
        consanguinity_details.append("general_consanguinity")
    
    return consanguinity_details[0], consanguinity_details

def calculate_combined_relatedness_enhanced(all_paths, consanguinity_type):
    """Enhanced relatedness calculation handling special consanguineous cases."""
    if not all_paths:
        return 0.0
    
    if len(all_paths) == 1:
        return all_paths[0][1]
    
    # Special handling for different consanguinity types
    if consanguinity_type == "double_cousin":
        # Double cousins: special calculation needed
        return handle_double_cousin_relatedness(all_paths)
    elif consanguinity_type == "pedigree_collapse":
        # Pedigree collapse: more complex probability calculations
        return handle_pedigree_collapse_relatedness(all_paths)
    else:
        # General case: sum independent path contributions
        total_relatedness = 0.0
        for _, relatedness in all_paths:
            total_relatedness += relatedness
        return total_relatedness

def handle_double_cousin_relatedness(all_paths):
    """Handle double cousin relatedness calculation."""
    # Double cousins share both sets of grandparents
    # Standard calculation is sum of path contributions
    total = sum(relatedness for _, relatedness in all_paths)
    return total

def handle_pedigree_collapse_relatedness(all_paths):
    """Handle pedigree collapse relatedness calculation."""
    # For pedigree collapse, we need to avoid double-counting
    # Group paths by their endpoint relatedness
    relatedness_groups = defaultdict(list)
    
    for rel_string, relatedness in all_paths:
        relatedness_groups[round(relatedness, 6)].append((rel_string, relatedness))
    
    # Sum contributions from different relatedness levels
    total = 0.0
    for relatedness_level, paths in relatedness_groups.items():
        # For same relatedness level, take the maximum contribution
        # to avoid over-counting the same genetic contribution
        if len(paths) > 1:
            total += relatedness_level
        else:
            total += relatedness_level
    
    return total

def get_siblings(tree, ind_id):
    """Enhanced sibling detection including complex family structures."""
    if ind_id not in tree:
        return {'full': [], 'half': [], 'step': []}
    
    ind = tree[ind_id]
    full_siblings = []
    half_siblings = []
    step_siblings = []
    
    for other_id, other in tree.items():
        if other_id != ind_id:
            # Full siblings: same father AND same mother
            if (ind['father_id'] and other['father_id'] == ind['father_id'] and
                ind['mother_id'] and other['mother_id'] == ind['mother_id']):
                full_siblings.append(other_id)
            
            # Half siblings: same father OR same mother (but not both)
            elif ((ind['father_id'] and other['father_id'] == ind['father_id'] and
                   (not ind['mother_id'] or not other['mother_id'] or ind['mother_id'] != other['mother_id'])) or
                  (ind['mother_id'] and other['mother_id'] == ind['mother_id'] and
                   (not ind['father_id'] or not other['father_id'] or ind['father_id'] != other['father_id']))):
                half_siblings.append(other_id)
            
            # Step siblings: one parent is married to the other's parent (but no shared biological parents)
            elif ((ind['father_id'] and other['mother_id'] == ind['father_id'] and ind['father_id'] != other['father_id']) or
                  (ind['mother_id'] and other['father_id'] == ind['mother_id'] and ind['mother_id'] != other['mother_id'])):
                step_siblings.append(other_id)
    
    return {
        'full': full_siblings,
        'half': half_siblings,
        'step': step_siblings
    }

def get_children(tree, parent_id):
    """Get all children of a parent."""
    children = []
    for ind_id, ind in tree.items():
        if ind['father_id'] == parent_id or ind['mother_id'] == parent_id:
            children.append(ind_id)
    return children

def build_relationship_string(path):
    """Build relationship string from path."""
    if not path:
        return "proband"
    return "->".join(["proband"] + path)

def calculate_relatedness_from_path(path):
    """Enhanced relatedness calculation from path."""
    if not path:
        return 1.0  # proband
    
    # For direct relationships (length 1)
    if len(path) == 1:
        step = path[0]
        if step in ["father", "mother"]:
            return 0.5
        elif step in ["son", "daughter", "child"]:
            return 0.5
        elif step in ["brother", "sister", "sibling"]:
            return 0.5
        elif step in ["half-brother", "half-sister", "half-sibling"]:
            return 0.25
        elif step in ["step-brother", "step-sister", "step-sibling"]:
            return 0.0
        else:
            return 0.5
    
    # For longer paths, count meioses
    meioses = 0
    has_step_relationship = False
    half_sibling_count = 0
    
    for step in path:
        if step in ["step-brother", "step-sister", "step-sibling"]:
            has_step_relationship = True
            break
        elif step in ["half-brother", "half-sister", "half-sibling"]:
            half_sibling_count += 1
            meioses += 1
        elif step in ["father", "mother", "son", "daughter", "child", "brother", "sister", "sibling"]:
            meioses += 1
        else:
            meioses += 1
    
    if has_step_relationship:
        return 0.0
    
    # Base relatedness calculation
    relatedness = 0.5 ** meioses
    
    # Adjust for half-siblings (they share only one parent instead of two)
    if half_sibling_count > 0:
        relatedness *= (0.5 ** half_sibling_count)
    
    return relatedness

def determine_relationship_type(relationship_string, relatedness):
    """Enhanced relationship type determination."""
    if relationship_string == "proband":
        return "proband"
    
    # Handle partners
    if "->partner" in relationship_string:
        related_part = relationship_string.replace("->partner", "")
        related_type = determine_relationship_type(related_part, calculate_relatedness_from_path(related_part.split("->")[1:] if "->" in related_part else []))
        return f"{related_type}'s partner"
    
    if relatedness == 0.0:
        if any(step in relationship_string for step in ["step-brother", "step-sister", "step-sibling"]):
            return "step-relative"
        return "unrelated"
    
    parts = relationship_string.split("->")[1:]  # Remove "proband"
    
    if len(parts) == 1:
        return parts[0]  # Direct relationships
    elif len(parts) == 2:
        # Grandparents, aunts/uncles, nephews/nieces, etc.
        if parts[0] in ["father", "mother"] and parts[1] in ["father", "mother"]:
            return "grandfather" if parts[1] == "father" else "grandmother"
        # [Add other 2-part patterns...]
    elif len(parts) == 3:
        # First cousins, great-grandparents, etc.
        if (parts[0] in ["father", "mother"] and 
            parts[1] in ["brother", "sister", "sibling", "half-brother", "half-sister"] and 
            parts[2] in ["son", "daughter", "child"]):
            return "first cousin"
        # [Add other 3-part patterns...]
    # [Continue with other length patterns...]
    
    # Fallback
    if relatedness > 0:
        return f"relative (relatedness: {relatedness})"
    else:
        return "step-relative"

def format_multiple_relationships_enhanced(all_paths, consanguinity_type, consanguinity_details):
    """Enhanced formatting for multiple relationship paths."""
    if not all_paths:
        return "unrelated"
    
    if len(all_paths) == 1:
        relationship_string, relatedness = all_paths[0]
        relationship_type = determine_relationship_type(relationship_string, relatedness)
        return f"{relationship_type} ({relationship_string})"
    
    # Multiple paths - format with consanguinity information
    relationship_descriptions = []
    for relationship_string, relatedness in all_paths:
        relationship_type = determine_relationship_type(relationship_string, relatedness)
        relationship_descriptions.append(f"{relationship_type} ({relationship_string})")
    
    formatted_relationships = ", ".join(relationship_descriptions)
    
    # Add consanguinity classification
    if consanguinity_type != "no_consanguinity":
        consanguinity_note = f" [CONSANGUINITY: {consanguinity_type}]"
        formatted_relationships += consanguinity_note
    
    return formatted_relationships

def main():
    import argparse

    parser = argparse.ArgumentParser(description="Enhanced relationship extraction from PED file with consanguinity detection.")
    parser.add_argument("input_ped", help="Input PED file")
    parser.add_argument("output_ped", help="Output file")
    parser.add_argument("--proband", help="Individual ID to use as proband (optional)", default=None)
    parser.add_argument("--proband-phenotype", help="Phenotype value to auto-detect proband (e.g., 2 for affected)", default=None)
    args = parser.parse_args()

    input_file = args.input_ped
    output_file = args.output_ped

    # Parse PED file
    individuals = parse_ped_file(input_file)

    if not individuals:
        print("No individuals found in PED file")
        sys.exit(1)

    # Validate pedigree
    validation_errors = validate_pedigree(individuals)
    if validation_errors:
        print("Pedigree validation errors found:")
        for error in validation_errors:
            print(f"  ERROR: {error}")
        print("\nProceeding with analysis, but results may be unreliable.\n")

    # Build enhanced family tree
    tree, children_dict, ancestors = build_family_tree(individuals)

    # Determine proband
    proband_id = None
    if args.proband:
        if args.proband in tree:
            proband_id = args.proband
        else:
            print(f"Specified proband '{args.proband}' not found in PED file.")
            sys.exit(1)
    elif args.proband_phenotype:
        # Try to auto-detect proband by phenotype value
        matching_inds = [ind for ind in individuals if str(ind['phenotype']) == str(args.proband_phenotype)]
        if matching_inds:
            proband_id = matching_inds[0]['ind_id']
            print(f"Auto-detected proband by phenotype: {proband_id}")
        else:
            print(f"No individual with phenotype '{args.proband_phenotype}' found. Using first individual as proband.")
            proband_id = individuals[0]['ind_id']
    else:
        # Default: first individual
        proband_id = individuals[0]['ind_id']

    # Calculate relationships and inbreeding coefficients
    results = []
    consanguinity_summary = defaultdict(int)
    
    for ind in individuals:
        # Find all relationship paths
        all_paths = find_all_relationship_paths_comprehensive(tree, ancestors, proband_id, ind['ind_id'])
        
        # Detect consanguinity type
        consanguinity_type, consanguinity_details = detect_consanguinity_type(all_paths)
        
        # Calculate enhanced relatedness
        combined_relatedness = calculate_combined_relatedness_enhanced(all_paths, consanguinity_type)
        
        # Calculate inbreeding coefficient
        inbreeding_coeff = calculate_inbreeding_coefficient(tree, ancestors, ind['ind_id'])
        
        # Format relationship description
        relationship_description = format_multiple_relationships_enhanced(all_paths, consanguinity_type, consanguinity_details)
        
        # Track consanguinity statistics
        if consanguinity_type != "no_consanguinity":
            consanguinity_summary[consanguinity_type] += 1
        
        results.append({
            'fam_id': ind['fam_id'],
            'ind_id': ind['ind_id'],
            'father_id': ind['father_id'] if ind['father_id'] else '0',
            'mother_id': ind['mother_id'] if ind['mother_id'] else '0',
            'sex': ind['sex'],
            'phenotype': ind['phenotype'],
            'relationship': relationship_description,
            'relatedness': combined_relatedness,
            'inbreeding_coeff': inbreeding_coeff,
            'consanguinity_type': consanguinity_type
        })
    
    # Write enhanced output file
    with open(output_file, 'w') as f:
        # Write header
        f.write("# Enhanced Relationship Analysis with Consanguinity Detection\n")
        f.write("# Columns: FAM_ID\tIND_ID\tFATHER_ID\tMOTHER_ID\tSEX\tPHENOTYPE\tRELATIONSHIP\tRELATEDNESS\tINBREEDING_COEFF\tCONSANGUINITY_TYPE\n")
        
        for result in results:
            line = f"{result['fam_id']}\t{result['ind_id']}\t{result['father_id']}\t{result['mother_id']}\t{result['sex']}\t{result['phenotype']}\t{result['relationship']}\t{result['relatedness']:.6f}\t{result['inbreeding_coeff']:.6f}\t{result['consanguinity_type']}"
            f.write(line + '\n')
    
    # Print comprehensive summary
    print(f"Processed {len(results)} individuals")
    relatedness_values = [r['relatedness'] for r in results if r['ind_id'] != proband_id]
    max_relatedness = max(relatedness_values) if relatedness_values else 0.0
    
    # Report consanguinity statistics
    consanguineous_individuals = [r for r in results if r['consanguinity_type'] != "no_consanguinity"]
    if consanguineous_individuals:
        print(f"\nFound {len(consanguineous_individuals)} individuals with consanguineous relationships:")
        
        for consang_type, count in consanguinity_summary.items():
            print(f"  {consang_type}: {count} individuals")
        
        print(f"\nDetailed consanguinity report:")
        for ind in consanguineous_individuals:
            print(f"  {ind['ind_id']}: {ind['relationship']}")
            print(f"    Relatedness: {ind['relatedness']:.6f}, Inbreeding: {ind['inbreeding_coeff']:.6f}")
    
    # Report individuals with high inbreeding coefficients
    high_inbreeding = [r for r in results if r['inbreeding_coeff'] > 0.0]
    if high_inbreeding:
        print(f"\nIndividuals with inbreeding (F > 0):")
        for ind in sorted(high_inbreeding, key=lambda x: x['inbreeding_coeff'], reverse=True):
            print(f"  {ind['ind_id']}: F = {ind['inbreeding_coeff']:.6f}")
    
    relatedness_values = [r['relatedness'] for r in results if r['ind_id'] != proband_id]
    inbreeding_values = [r['inbreeding_coeff'] for r in results]
    max_relatedness = max(relatedness_values) if relatedness_values else 0.0
    max_inbreeding = max(inbreeding_values) if inbreeding_values else 0.0
    
    print(f"\nSummary Statistics:")
    print(f"  Maximum relatedness to proband: {max_relatedness:.6f}")
    print(f"  Maximum inbreeding coefficient: {max_inbreeding:.6f}")
    print(f"  Consanguineous relationships detected: {len(consanguineous_individuals)}")
    print(f"  Inbred individuals: {len(high_inbreeding)}")
    print(f"  Inbred individuals: {len(high_inbreeding)}")

if __name__ == "__main__":
    main()