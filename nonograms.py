import matplotlib.pyplot as plt
import random
import numpy as np
import json
import math
import pygad
import pyswarms as ps


class Nonogram:
    def __init__(self, columns, rows, sol=None):
        self.columns = columns
        self.rows = rows
        self.sol = sol if sol is not None else self.generate_empty_solution()

    def __str__(self):
        return f"Columns: {self.columns}\nRows: {self.rows}\nSolution: {self.sol}"

    def verify(self):
        def check_line(line, pattern):
            count = []
            current = 0
            for cell in line:
                if cell == 1:
                    current += 1
                elif current > 0:
                    count.append(current)
                    current = 0
            if current > 0:
                count.append(current)
            return count == pattern

        for i, row in enumerate(self.sol):
            if not check_line(row, self.rows[i]):
                return False

        for j in range(len(self.columns)):
            col = [self.sol[i][j] for i in range(len(self.sol))]
            if not check_line(col, self.columns[j]):
                return False

        return True

    def visualize(self):
        fig, ax = plt.subplots()
        ax.set_xticks([])
        ax.set_yticks([])

        for i in range(len(self.sol) + 1):
            ax.plot([0, len(self.sol)], [i, i], "k-")
        for j in range(len(self.sol[0]) + 1):
            ax.plot([j, j], [0, len(self.sol)], "k-")

        for i in range(len(self.sol)):
            for j in range(len(self.sol[0])):
                if self.sol[i][j] == 1:
                    ax.fill_between(
                        [j, j + 1],
                        len(self.sol) - i - 1,
                        len(self.sol) - i,
                        color="black",
                    )

        for j in range(len(self.columns)):
            col = self.columns[j]
            col_str = ", ".join(map(str, col))
            ax.text(
                j + 0.5,
                len(self.sol) + 0.5,
                col_str,
                color="black",
                ha="center",
                va="bottom",
            )

        for i in range(len(self.rows)):
            row = self.rows[i]
            row_str = ", ".join(map(str, row))
            ax.text(
                -0.5,
                len(self.sol) - i - 0.5,
                row_str,
                color="black",
                ha="right",
                va="center",
            )
        plt.grid(color="red", linestyle="-", linewidth=2)
        plt.show()

    def calculate_rows_columns(self, sol):
        rows = [
            [
                sum(group)
                for group in np.split(row, np.where(np.diff(row) != 0)[0] + 1)
                if group[0] != 0
            ]
            for row in sol
        ]
        columns = [
            [
                sum(group)
                for group in np.split(column, np.where(np.diff(column) != 0)[0] + 1)
                if group[0] != 0
            ]
            for column in zip(*sol)
        ]
        return rows, columns

    def generate_empty_solution(self):
        return [[0 for _ in range(len(self.columns))] for _ in range(len(self.rows))]

    def from_vector(self, vector):
        return np.array(vector).reshape(len(self.rows), len(self.columns))

    @staticmethod
    def generate_random_nonogram(n, m):
        sol = [[random.randint(0, 1) for _ in range(m)] for _ in range(n)]

        rows = [
            [
                sum(group)
                for group in np.split(row, np.where(np.diff(row) != 0)[0] + 1)
                if group[0] != 0
            ]
            for row in sol
        ]
        columns = [
            [
                sum(group)
                for group in np.split(column, np.where(np.diff(column) != 0)[0] + 1)
                if group[0] != 0
            ]
            for column in zip(*sol)
        ]

        nonogram = Nonogram(columns, rows, sol)

        calculated_rows, calculated_columns = nonogram.calculate_rows_columns(
            nonogram.sol
        )
        assert calculated_rows == nonogram.rows
        assert calculated_columns == nonogram.columns

        return nonogram

    def runGA(
        self,
        fitness_func,
        generations,
        sol_per_pop,
        mutation_percent_genes,
        num_parents_mating,
    ):
        @staticmethod
        def arrayEquals(arr1, arr2):
            if len(arr1) != len(arr2):
                return False
            for i in range(len(arr1)):
                if arr1[i] != arr2[i]:
                    return False
            return True

        def on_generation(ga_instance):
            best_solution, solution_fitness, solution_idx = ga_instance.best_solution()
            print(
                f"Generation {ga_instance.generations_completed}: Best Fitness = {solution_fitness}"
            )

        def fitness(ga_instance, solution, solution_idx):
            resultRow = []
            resultCol = []
            solution = [
                solution[i : i + len(self.rows)]
                for i in range(0, len(solution), len(self.rows))
            ]
            for i in solution:
                count = 0
                tmp = []
                for j in i:
                    if j == 1:
                        count += 1
                    elif count != 0:
                        tmp.append(count)
                        count = 0
                if count != 0:
                    tmp.append(count)
                resultRow.append(tmp)
            for i in range(len(solution[0])):
                count = 0
                tmp = []
                for j in range(len(solution)):
                    if solution[j][i] == 1:
                        count += 1
                    elif count != 0:
                        tmp.append(count)
                        count = 0
                if count != 0:
                    tmp.append(count)
                resultCol.append(tmp)
            wrong = 0
            for i in range(len(self.rows)):
                if not arrayEquals(self.rows[i], resultRow[i]):
                    wrong -= 1
                if not arrayEquals(self.columns[i], resultCol[i]):
                    wrong -= 1
            return wrong

        def fitness_advanced(ga_instance, solution, solution_idx):
            resultRow = []
            resultCol = []
            solution = [
                solution[i : i + len(self.rows)]
                for i in range(0, len(solution), len(self.rows))
            ]

            # Calculate row blocks
            for i in solution:
                count = 0
                tmp = []
                for j in i:
                    if j == 1:
                        count += 1
                    else:
                        if count != 0:
                            tmp.append(count)
                            count = 0
                if count != 0:
                    tmp.append(count)
                resultRow.append(tmp)

            # Calculate column blocks
            for i in range(len(solution[0])):
                count = 0
                tmp = []
                for j in range(len(solution)):
                    if solution[j][i] == 1:
                        count += 1
                    else:
                        if count != 0:
                            tmp.append(count)
                            count = 0
                if count != 0:
                    tmp.append(count)
                resultCol.append(tmp)

            wrong = 0

            # Evaluate rows
            for i in range(len(self.rows)):
                if i >= len(resultRow):
                    wrong -= 10 * abs(len(self.rows[i]))  # Penalty for missing rows
                    continue
                if len(self.rows[i]) != len(resultRow[i]):
                    wrong -= (
                        abs(len(self.rows[i]) - len(resultRow[i])) * 10
                    )  # Larger penalty for block count difference
                min_len = min(len(self.rows[i]), len(resultRow[i]))
                for j in range(min_len):
                    wrong -= abs(
                        self.rows[i][j] - resultRow[i][j]
                    )  # Penalty for length mismatch
                # Extra penalty for position mismatch
                for j in range(min_len):
                    if self.rows[i][j] != resultRow[i][j]:
                        wrong -= 1

            # Evaluate columns
            for i in range(len(self.columns)):
                if i >= len(resultCol):
                    wrong -= 10 * abs(
                        len(self.columns[i])
                    )  # Penalty for missing columns
                    continue
                if len(self.columns[i]) != len(resultCol[i]):
                    wrong -= (
                        abs(len(self.columns[i]) - len(resultCol[i])) * 10
                    )  # Larger penalty for block count difference
                min_len = min(len(self.columns[i]), len(resultCol[i]))
                for j in range(min_len):
                    wrong -= abs(
                        self.columns[i][j] - resultCol[i][j]
                    )  # Penalty for length mismatch
                # Extra penalty for position mismatch
                for j in range(min_len):
                    if self.columns[i][j] != resultCol[i][j]:
                        wrong -= 1

            return wrong

        def fitness_another(ga_instance, solution, solution_idx):
            tab = [[0 for _ in range(len(self.columns))] for _ in range(len(self.rows))]
            iterator = 0
            wrong = 0
            solution = solution.astype(int)

            # Process rows
            for i in range(len(self.rows)):
                for j in range(len(self.rows[i])):
                    for k in range(self.rows[i][j]):
                        if iterator + k < len(solution):
                            if solution[iterator + k] < len(self.rows):
                                tab[i][solution[iterator + k]] = 1
                            else:
                                wrong -= 1
                        else:
                            wrong -= 1
                    iterator += 1

            # Process columns
            for i in range(len(self.columns)):
                for j in range(len(self.columns[i])):
                    for k in range(self.columns[i][j]):
                        if iterator + k < len(solution):
                            if solution[iterator + k] < len(self.columns):
                                tab[solution[iterator + k]][i] = 1
                            else:
                                wrong -= 1
                        else:
                            wrong -= 1
                    iterator += 1

            resultRow = []
            resultCol = []

            # Process result rows
            for i in tab:
                counter = 0
                temp = []
                for j in i:
                    if j == 1:
                        counter += 1
                    elif counter != 0:
                        temp.append(counter)
                        counter = 0
                if counter != 0:
                    temp.append(counter)
                resultRow.append(temp)

            # Process result columns
            for i in range(len(tab[0])):
                counter = 0
                temp = []
                for j in range(len(tab)):
                    if tab[j][i] == 1:
                        counter += 1
                    elif counter != 0:
                        temp.append(counter)
                        counter = 0
                if counter != 0:
                    temp.append(counter)
                resultCol.append(temp)

            # Evaluate rows
            for i in range(len(self.rows)):
                if i < len(resultRow):
                    if len(self.rows[i]) == len(resultRow[i]):
                        for j in range(len(self.rows[i])):
                            if self.rows[i][j] != resultRow[i][j]:
                                wrong -= 1
                    else:
                        wrong -= abs(len(self.rows[i]) - len(resultRow[i]))
                        if len(self.rows[i]) > len(resultRow[i]):
                            for j in range(len(resultRow[i])):
                                if self.rows[i][j] != resultRow[i][j]:
                                    wrong -= 1
                        else:
                            for j in range(len(self.rows[i])):
                                if self.rows[i][j] != resultRow[i][j]:
                                    wrong -= 1
                else:
                    wrong -= 10  # Penalty for missing rows

            # Evaluate columns
            for i in range(len(self.columns)):
                if i < len(resultCol):
                    if len(self.columns[i]) == len(resultCol[i]):
                        for j in range(len(self.columns[i])):
                            if self.columns[i][j] != resultCol[i][j]:
                                wrong -= 1
                    else:
                        wrong -= abs(len(self.columns[i]) - len(resultCol[i]))
                        if len(self.columns[i]) > len(resultCol[i]):
                            for j in range(len(resultCol[i])):
                                if self.columns[i][j] != resultCol[i][j]:
                                    wrong -= 1
                        else:
                            for j in range(len(self.columns[i])):
                                if self.columns[i][j] != resultCol[i][j]:
                                    wrong -= 1
                else:
                    wrong -= 10  # Penalty for missing columns

            return wrong

        def best_fit(ga_instance, solution, solution_idx):
            tab = [[0] * len(self.columns) for _ in range(len(self.rows))]
            iterator = 0
            wrong = 0
            solution = solution.astype(int)

            for i, row in enumerate(self.rows):
                for j, count in enumerate(row):
                    for k in range(count):
                        pos = solution[iterator] + k
                        if pos < len(self.rows):
                            tab[i][pos] = 1
                        else:
                            wrong -= 1
                    iterator += 1
            resultCol = []
            for col in zip(*tab):
                temp = []
                counter = 0
                for cell in col:
                    if cell == 1:
                        counter += 1
                    elif counter:
                        temp.append(counter)
                        counter = 0
                if counter:
                    temp.append(counter)
                resultCol.append(temp)

            for col, res_col in zip(self.columns, resultCol):
                if len(col) == len(res_col):
                    wrong -= sum(1 for a, b in zip(col, res_col) if a != b)
                else:
                    wrong -= abs(len(col) - len(res_col))
                    if len(col) > len(res_col):
                        wrong -= sum(1 for a, b in zip(col, res_col) if a != b)
                    else:
                        wrong -= sum(1 for a, b in zip(res_col, col) if a != b)

            return wrong

        gene_range = [0, 1]
        num_genes = len(self.rows) * len(self.columns)

        if fitness_func == "fitness":
            fit = fitness
        elif fitness_func == "fitness_advanced":
            fit = fitness_advanced
        elif fitness_func == "fitness_another":
            fit = fitness_another
            gene_range = [i for i in range(len(self.rows))]
            num_genes = 0
            for i in range(len(self.rows)):
                num_genes += len(self.rows[i])
                num_genes += len(self.columns[i])
        elif fitness_func == "best_fit":
            fit = best_fit
            gene_range = [i for i in range(len(self.rows))]
            num_genes = 0
            for i in range(len(self.rows)):
                num_genes += len(self.rows[i])

        num_generations = generations
        keep_parents = 2
        parent_type_selection = "sss"
        crossover_type = "single_point"
        mutation_type = "random"

        ga_instance = pygad.GA(
            gene_space=gene_range,
            num_generations=num_generations,
            num_parents_mating=num_parents_mating,
            fitness_func=fit,
            sol_per_pop=sol_per_pop,
            num_genes=num_genes,
            parent_selection_type=parent_type_selection,
            keep_parents=keep_parents,
            crossover_type=crossover_type,
            mutation_type=mutation_type,
            mutation_percent_genes=mutation_percent_genes,
            on_generation=on_generation,
        )

        ga_instance.run()

        solution, solution_fitness, solution_idx = ga_instance.best_solution()

        ga_instance.plot_fitness()
        ga_instance.save("ga_fitness")

        return {
            "best_solution_generation": ga_instance.best_solution_generation,
            "solution_fitness": solution_fitness,
            "solution": solution,
        }

    def runSwarm(self, fitness_func, generations, c1, c2, w, k=None, p=None):
        def arrayEquals(arr1, arr2):
            if len(arr1) != len(arr2):
                return False
            for i in range(len(arr1)):
                if arr1[i] != arr2[i]:
                    return False
            return True

        def fitness(temp):
            resultRow = []
            resultCol = []
            solution = [i[0] for i in temp]
            solution = [
                solution[i : i + len(self.rows)]
                for i in range(0, len(solution), len(self.rows))
            ]

            for i in solution:
                counter = 0
                temp = []
                for j in i:
                    if j == 1:
                        counter += 1
                    elif counter != 0:
                        temp.append(counter)
                        counter = 0
                if counter != 0:
                    temp.append(counter)
                resultRow.append(temp)

            for i in range(len(solution[0])):
                counter = 0
                temp = []
                for j in range(len(solution[i])):
                    if solution[j][i] == 1:
                        counter += 1
                    elif counter != 0:
                        temp.append(counter)
                        counter = 0
                if counter != 0:
                    temp.append(counter)
                resultCol.append(temp)
            wrongCounter = 0
            for i in range(len(self.rows)):
                if len(self.rows[i]) == len(resultRow[i]):
                    for j in range(len(self.rows[i])):
                        if self.rows[i][j] != resultRow[i][j]:
                            wrongCounter -= 1
                else:
                    wrongCounter -= math.fabs(len(self.rows[i]) - len(resultRow[i]))
                    if len(self.rows[i]) > len(resultRow[i]):
                        for j in range(len(resultRow[i])):
                            if self.rows[i][j] != resultRow[i][j]:
                                wrongCounter -= 1
                    else:
                        for j in range(len(self.rows[i])):
                            if self.rows[i][j] != resultRow[i][j]:
                                wrongCounter -= 1
                if len(self.columns[i]) == len(resultCol[i]):
                    for j in range(len(self.columns[i])):
                        if self.columns[i][j] != resultCol[i][j]:
                            wrongCounter -= 1
                else:
                    wrongCounter -= math.fabs(len(self.columns[i]) - len(resultCol[i]))
                    if len(self.columns[i]) > len(resultCol[i]):
                        for j in range(len(resultCol[i])):
                            if self.columns[i][j] != resultCol[i][j]:
                                wrongCounter -= 1
                    else:
                        for j in range(len(self.columns[i])):
                            if self.columns[i][j] != resultCol[i][j]:
                                wrongCounter -= 1
            return -wrongCounter

        def best_fit(solution):
            # Inicjalizacja siatki
            tab = np.zeros((len(self.rows), len(self.columns)), dtype=int)
            wrong_counter = 0

            # Konwersja solution na listę liczb całkowitych
            solution = np.round(solution).astype(int).flatten().tolist()

            # Wypełnianie siatki na podstawie solution
            iterator = 0
            for row_idx, row in enumerate(self.rows):
                for segment_length in row:
                    start_position = solution[iterator]
                    if start_position + segment_length <= len(self.columns):
                        tab[
                            row_idx, start_position : start_position + segment_length
                        ] = 1
                    else:
                        wrong_counter -= 1
                    iterator += 1

            # Generowanie wynikowych kolumn
            result_col = []
            for col_idx in range(tab.shape[1]):
                col_segment_lengths = []
                counter = 0
                for row_idx in range(tab.shape[0]):
                    if tab[row_idx, col_idx] == 1:
                        counter += 1
                    elif counter > 0:
                        col_segment_lengths.append(counter)
                        counter = 0
                if counter > 0:
                    col_segment_lengths.append(counter)
                result_col.append(col_segment_lengths)

            # Sprawdzanie poprawności kolumn
            for col_idx, col in enumerate(self.columns):
                if len(col) != len(result_col[col_idx]):
                    wrong_counter -= abs(len(col) - len(result_col[col_idx]))
                for seg_idx, seg_length in enumerate(col):
                    if (
                        seg_idx >= len(result_col[col_idx])
                        or seg_length != result_col[col_idx][seg_idx]
                    ):
                        wrong_counter -= 1

            return -wrong_counter

        options = {"c1": c1, "c2": c2, "w": w, "k": k, "p": p}
        if fitness_func == "fitness":
            fitness_func = fitness

            n = len(self.rows) * len(self.columns)
            optimizer = ps.discrete.BinaryPSO(
                n_particles=30, dimensions=n, options=options
            )

            final_sol = {
                "best_solution_generation": None,
                "solution_fitness": None,
                "solution": None,
            }

            stats = optimizer.optimize(fitness_func, iters=generations, verbose=False)

            cost_history = optimizer.cost_history

            final_sol["solution_fitness"] = optimizer.cost_history[
                len(optimizer.cost_history) - 1
            ]
            best_solution = optimizer.swarm.position[
                np.argmin(optimizer.swarm.best_cost)
            ]
            final_sol["solution"] = np.round(best_solution).reshape(
                len(self.rows), len(self.columns)
            )
            for i in range(len(cost_history)):
                if cost_history[i] == final_sol["solution_fitness"]:
                    final_sol["best_solution_generation"] = i
                    break
            return final_sol
        else:
            limiter = None
            num_genes = None
            fitness_func = best_fit
            limiter = len(self.rows)
            num_genes = 0
            for i in range(len(self.rows)):
                num_genes += len(self.rows[i])

            options = {"c1": 0.5, "c2": 0.3, "w": 0.9}

            max_bound = limiter * np.ones(num_genes)
            min_bound = 0 * np.ones(num_genes)
            bounds = (min_bound, max_bound)

            optimizer = ps.single.GlobalBestPSO(
                n_particles=30, dimensions=num_genes, options=options, bounds=bounds
            )

            final_sol = {
                "best_solution_generation": None,
                "solution_fitness": None,
                "solution": None,
            }

            stats = optimizer.optimize(fitness_func, iters=generations, verbose=False)

            cost_history = optimizer.cost_history

            final_sol["solution_fitness"] = optimizer.cost_history[
                len(optimizer.cost_history) - 1
            ]
            best_solution = optimizer.swarm.position[
                np.argmin(optimizer.swarm.best_cost)
            ]
            final_sol["solution"] = best_solution
            for i in range(len(cost_history)):
                if cost_history[i] == final_sol["solution_fitness"]:
                    final_sol["best_solution_generation"] = i
                    break
            return final_sol

    @staticmethod
    def vector_to_binary_grid(solution, rows, columns):
        grid = [[0 for _ in range(len(columns))] for _ in range(len(rows))]

        iterator = 0
        for i in range(len(rows)):
            for j in range(len(rows[i])):
                start_position = int(solution[iterator])
                length = rows[i][j]
                for k in range(length):
                    if start_position + k < len(grid[i]):
                        grid[i][start_position + k] = 1
                iterator += 1
        print(grid)
        return np.array(grid).flatten()
