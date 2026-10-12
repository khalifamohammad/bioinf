import random
from collections import Counter


def generate_sequence(length=20):
    """Generate a random DNA sequence."""
    return ''.join(random.choice("ACTG") for _ in range(length))


def nucleotide_frequency(sequence):
    """Find frequency of individual nucleotides."""
    print("\n Single Letter Frequencies ")

    counts = Counter(sequence)

    for base in "ACTG":
        count = counts[base]
        frequency = count / len(sequence)

        print(f"{base}: {count} ({frequency:.2%})")


def combination_frequency(sequence, length):
    """Find frequencies of overlapping combinations."""

    combinations = []

    for i in range(len(sequence) - length + 1):
        combination = sequence[i:i + length]
        combinations.append(combination)

    counts = Counter(combinations)

    sorted_counts = counts.most_common()

    print(f"\n {length}-Letter Combinations ")

    total = len(combinations)

    for combination, count in sorted_counts:
        frequency = count / total

        print(f"{combination}: {count} ({frequency:.2%})")




sequence = generate_sequence(20)

print("Random DNA Sequence:")
print(sequence)

nucleotide_frequency(sequence)


combination_frequency(sequence, 2)


combination_frequency(sequence, 3)