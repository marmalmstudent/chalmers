def func(old_val, orig_pwr, x, n):
    if (n > 1):
        return func(old_val + orig_pwr*x**n, orig_pwr, x, n-1)
    elif (n == 0):
        return old_val
    else:
        return old_val + orig_pwr*x**n


n = 280
print("%.1f" % func(200, 200, 0.99, n), " W after %.0f" % n, " bounces")
