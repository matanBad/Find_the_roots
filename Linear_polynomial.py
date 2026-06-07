def linear_interpolation(x_points, y_points, x_target):
    if x_target in x_points:
        return y_points[x_points.index(x_target)]

    if x_target < min(x_points) or x_target > max(x_points):
        print("The target x is outside the table range.")
        return None

    for i in range(len(x_points) - 1):
        x1 = x_points[i]
        x2 = x_points[i + 1]
        y1 = y_points[i]
        y2 = y_points[i + 1]

        if x1 <= x_target <= x2:
            return y1 + ((y2 - y1) / (x2 - x1)) * (x_target - x1)

    return None


def solve_linear_system(matrix, vector):
    n = len(matrix)

    for i in range(n):
        max_row = i
        for k in range(i + 1, n):
            if abs(matrix[k][i]) > abs(matrix[max_row][i]):
                max_row = k

        matrix[i], matrix[max_row] = matrix[max_row], matrix[i]
        vector[i], vector[max_row] = vector[max_row], vector[i]

        if matrix[i][i] == 0:
            print("Cannot solve system: zero pivot.")
            return None

        for k in range(i + 1, n):
            factor = matrix[k][i] / matrix[i][i]

            for j in range(i, n):
                matrix[k][j] -= factor * matrix[i][j]

            vector[k] -= factor * vector[i]

    solution = [0.0] * n

    for i in range(n - 1, -1, -1):
        total = 0.0

        for j in range(i + 1, n):
            total += matrix[i][j] * solution[j]

        solution[i] = (vector[i] - total) / matrix[i][i]

    return solution


def polynomial_interpolation(x_points, y_points, x_target):
    if x_target in x_points:
        return y_points[x_points.index(x_target)]

    if x_target < min(x_points) or x_target > max(x_points):
        print("The target x is outside the table range.")
        return None

    n = len(x_points)
    matrix = []

    for x in x_points:
        row = []
        for power in range(n):
            row.append(x ** power)
        matrix.append(row)

    vector = y_points[:]

    coefficients = solve_linear_system(matrix, vector)

    if coefficients is None:
        return None

    result = 0.0

    for power in range(n):
        result += coefficients[power] * (x_target ** power)

    return result


def main():
    x_points = [1.0, 2.0, 3.0, 4.0]
    y_points = [1.0, 4.0, 9.0, 16.0]

    x_target = 2.5

    print("Interpolation Program")
    print("---------------------")

    print("Table points:")
    for i in range(len(x_points)):
        print(f"({x_points[i]}, {y_points[i]})")

    print(f"\nTarget x value: {x_target}")

    linear_result = linear_interpolation(x_points, y_points, x_target)
    polynomial_result = polynomial_interpolation(x_points, y_points, x_target)

    print("\nResults:")

    if linear_result is not None:
        print(f"Linear interpolation result: y ≈ {linear_result}")

    if polynomial_result is not None:
        print(f"Polynomial interpolation result: y ≈ {polynomial_result}")


if __name__ == "__main__":
    main()