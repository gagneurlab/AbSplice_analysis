from scipy.stats import fisher_exact


def contingency_table(first_set: set, second_set: set, universal_set: set):
    return [
        [
            len(first_set.intersection(second_set)),
            len(first_set.difference(second_set))
        ],
        [
            len(second_set.difference(first_set)),
            len(universal_set.difference(first_set.union(second_set)))
        ]
    ]


def fisher_exact_sets(first_set: set, second_set: set, universal_set: set):
    return fisher_exact(contingency_table(first_set, second_set, universal_set))
