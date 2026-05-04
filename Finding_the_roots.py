import sympy


def initializeSympyPolynomialData(poly_expr, symbol):
    print("Function:", poly_expr)

    # חישוב הנגזרת הראשונה והשנייה
    derivative_expr = sympy.diff(poly_expr, symbol)
    second_derivative_expr = sympy.diff(derivative_expr, symbol)

    f = sympy.utilities.lambdify(symbol, poly_expr, 'math')
    fTag = sympy.utilities.lambdify(symbol, derivative_expr, 'math')
    fTagTag = sympy.utilities.lambdify(symbol, second_derivative_expr, 'math')

    return f, fTag, fTagTag


def Bisection_Method(poly_func, start_point, end_point, epsilon=0.0001):
    if poly_func(start_point) * poly_func(end_point) >= 0:
        return None, 0

    print(f"\n--- Searching in range [{start_point:.2f}, {end_point:.2f}] ---")
    print("| Iteration | Current x (Root) |")
    print("|-----------|------------------|")

    iterations = 0
    a, b = start_point, end_point
    c = a

    while (b - a) >= epsilon:
        iterations += 1
        c = (a + b) / 2.0
        print(f"| {iterations:<9} | {c:<16.5f} |")

        if poly_func(c) == 0.0:
            break
        elif poly_func(c) * poly_func(a) < 0:
            b = c
        else:
            a = c

    print("-" * 32)
    print(f"-> Final Root found: {c:.5f} | Total Iterations for this root: {iterations}")
    return c, iterations


def Newton_Raphson(poly_func, derivative_func, start_point, end_point, epsilon=0.0001):
    current_x = (start_point + end_point) / 2.0
    iterations = 0
    max_iterations = 100

    print(f"\n--- Searching in range [{start_point:.2f}, {end_point:.2f}] ---")
    print("| Iteration | Current x (Root) |")
    print("|-----------|------------------|")

    while iterations < max_iterations:
        f_x = poly_func(current_x)
        f_prime_x = derivative_func(current_x)

        if f_prime_x == 0: return None, 0

        next_x = current_x - (f_x / f_prime_x)
        iterations += 1
        print(f"| {iterations:<9} | {next_x:<16.5f} |")

        if abs(next_x - current_x) < epsilon:
            if min(start_point, end_point) <= next_x <= max(start_point, end_point):
                print("-" * 32)
                print(f"-> Final Root found: {next_x:.5f} | Total Iterations for this root: {iterations}")
                return next_x, iterations
            return None, 0

        current_x = next_x

    return None, 0


def secant_method(poly_func, start_point, end_point, epsilon=0.0001):
    x0, x1 = start_point, end_point
    iterations = 0
    max_iterations = 100

    print(f"\n--- Searching in range [{start_point:.2f}, {end_point:.2f}] ---")
    print("| Iteration | Current x (Root) |")
    print("|-----------|------------------|")

    while iterations < max_iterations:
        f_x0 = poly_func(x0)
        f_x1 = poly_func(x1)

        if abs(f_x1 - f_x0) < 1e-12: return None, 0

        x2 = x1 - f_x1 * ((x1 - x0) / (f_x1 - f_x0))
        iterations += 1
        print(f"| {iterations:<9} | {x2:<16.5f} |")

        if abs(x2 - x1) < epsilon or abs(poly_func(x2)) < epsilon:
            if min(start_point, end_point) - epsilon <= x2 <= max(start_point, end_point) + epsilon:
                print("-" * 32)
                print(f"-> Final Root found: {x2:.5f} | Total Iterations for this root: {iterations}")
                return x2, iterations
            return None, 0

        # התיקון: מעדכנים את שתי הנקודות לקראת האיטרציה הבאה
        x0, x1 = x1, x2

    return None, 0


# =========================================================
# התוכנית הראשית
# =========================================================
if __name__ == "__main__":
    x = sympy.symbols('x')
    my_polynomial = sympy.cos((2 * x ** 3) + (5 * x ** 2) - 6) / (2 * sympy.exp(-2 * x))

    print("Initializing Sympy Data...")
    f_func, fTag_func, fTagTag_func = initializeSympyPolynomialData(my_polynomial, x)

    # טווח הבדיקה והחלוקה למקטעים של 0.1
    large_start = -2.0
    large_end = 2.0
    step = 0.1

    print("\nSelect the method to find the roots:")
    print("1. Bisection Method")
    print("2. Newton-Raphson")
    print("3. Secant Method")
    choice = input("Enter your choice (1/2/3): ")

    found_roots = []
    current = large_start

    print("\nScanning for roots...")
    while current < large_end:
        a = current
        b = current + step
        if b > large_end: b = large_end

        # 1. בדיקה לחיתוך רגיל של ציר ה-x
        if f_func(a) * f_func(b) < 0:
            root, iters = None, 0
            if choice == '1':
                root, iters = Bisection_Method(f_func, a, b)
            elif choice == '2':
                root, iters = Newton_Raphson(f_func, fTag_func, a, b)
            elif choice == '3':
                root, iters = secant_method(f_func, a, b)

            if root is not None and not any(abs(root - r[0]) < 0.01 for r in found_roots):
                found_roots.append((root, iters, "Odd (Regular cross)"))

        # 2. בדיקה להשקה בציר ה-x (שורש מריבוב זוגי)
        if fTag_func(a) * fTag_func(b) < 0:
            ext_root, iters = None, 0
            if choice == '1':
                ext_root, iters = Bisection_Method(fTag_func, a, b)
            elif choice == '2':
                ext_root, iters = Newton_Raphson(fTag_func, fTagTag_func, a, b)
            elif choice == '3':
                ext_root, iters = secant_method(fTag_func, a, b)

            if ext_root is not None and abs(f_func(ext_root)) < 0.001:
                if not any(abs(ext_root - r[0]) < 0.01 for r in found_roots):
                    found_roots.append((ext_root, iters, "Even (Tangent touch)"))

        current += step

    # הדפסת טבלת התוצאות הסופית והנקייה
    print("\n--- Final Results ---")
    if not found_roots:
        print("No roots found in the given range.")
    else:
        print("| Root # | x          | Type                  |")
        print("|--------|------------|-----------------------|")
        for i, (r, iters, r_type) in enumerate(found_roots):
            print(f"| {i + 1:<6} | {r:<10.5f} | {r_type:<21} |")