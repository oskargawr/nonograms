import numpy as np

arr = [3, 3, 6, 0, 3, 0, 0, 3, 6, 3, 2]


def vector_to_binary_grid(solution, rows, columns):
    grid = [[0 for _ in range(len(columns))] for _ in range(len(rows))]

    print("grid:", grid)
    iterator = 0
    for i in range(len(rows)):
        for j in range(len(rows[i])):
            start_position = int(solution[iterator])
            length = rows[i][j]
            for k in range(length):
                if start_position + k < len(grid[i]):
                    grid[i][start_position + k] = 1

            iterator += 1
            # print("grid:", grid)
    print(grid)
    return np.array(grid).flatten()


print(
    vector_to_binary_grid(
        arr,
        [[3], [3], [1, 1], [7], [1, 3], [1, 1], [2, 2]],
        [[2], [1, 1], [2, 4], [5], [2, 4], [1, 1], [2]],
    )
)
