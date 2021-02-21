import numpy as np
import pytest

np.random.seed(665456)


@pytest.mark.parametrize("seed", np.random.randint(0, np.iinfo(np.uint32).max, 25))
@pytest.mark.parametrize("nsam", [i for i in range(10, 51)])
def test_update_R(seed, nsam):
    state = np.random.RandomState(seed=seed)
    R = np.arange(nsam)
    i = nsam
    X = 2 * nsam

    coalesced = []

    while i > 1:
        c0 = state.randint(0, i, 1)[0]
        c1 = state.randint(0, i, 1)[0]

        while c1 == c0:
            c1 = state.randint(0, i, 1)[0]

        coalesced.append(R[c0])
        coalesced.append(R[c1])

        R[min(c0, c1)] = X - i
        R[max(c0, c1)] = R[i - 1]

        i -= 1
    coalcounts = np.unique(coalesced, return_counts=True)
    for i in range(2 * nsam - 2):  # Skip the root
        assert coalcounts[1][i] == 1
