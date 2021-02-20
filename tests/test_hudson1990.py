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
    expected = {i: True for i in R}

    while i > 1:
        for r in R[:i]:
            assert r in expected
        c0 = state.randint(0, i, 1)[0]
        c1 = state.randint(0, i, 1)[0]

        while c1 == c0:
            c1 = state.randint(0, i, 1)[0]

        assert R[c0] in expected
        assert R[c1] in expected

        expected.pop(min(R[c0], R[c1]))
        expected.pop(max(R[c0], R[c1]))
        expected[X - i] = True
        expected[max(R[c0], R[c1])] = True

        R[min(c0, c1)] = X - i
        R[max(c0, c1)] = R[i - 1]

        i -= 1
