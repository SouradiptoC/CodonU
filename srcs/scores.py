from Bio.Data.IUPACData import unambiguous_dna_letters, unambiguous_rna_letters


def base_score(rna=False) -> dict:
    return {'U': 1, 'C': 2, 'A': 3, 'G': 4} if rna else {'T': 1, 'C': 2, 'A': 3, 'G': 4}


def codon_score(rna=False) -> dict:
    score_dict = dict()
    base_scores = base_score(rna)
    if rna:
        for i in unambiguous_rna_letters:
            for j in unambiguous_rna_letters:
                for k in unambiguous_rna_letters:
                    score_dict[''.join([i, j, k])] = ((base_scores[i] - 1) * 16) + base_scores[j] + (
                            (base_scores[k] - 1) * 4)
    else:
        for i in unambiguous_dna_letters:
            for j in unambiguous_dna_letters:
                for k in unambiguous_dna_letters:
                    score_dict[''.join([i, j, k])] = ((base_scores[i] - 1) * 16) + base_scores[j] + (
                            (base_scores[k] - 1) * 4)
    return score_dict


if __name__ == '__main__':
    scores = codon_score()
    print(scores['GCG'])
