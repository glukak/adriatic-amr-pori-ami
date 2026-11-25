"""Compares key numbers of the run with those reported in the manuscript."""
import sys
import pandas as pd

F = "results/AMR_results.xlsx"
N = pd.read_excel(F, sheet_name=None)
P = N["S3_PORI"].set_index("Location")
ami = N["S6s_AMI_source"]
sent = N["S4_sentinels_E"]
ok_all = True

def check(name, got, exp, tol):
    global ok_all
    ok = abs(got - exp) <= tol
    ok_all &= ok
    print(f"{'OK  ' if ok else 'FAIL'} {name}: got {got:.3f}, expected {exp} (tol {tol})")

tiers = P["PORI tier"].value_counts()
for t, n in [("Low", 22), ("Moderate", 8), ("High", 4), ("Very high", 1), ("Extreme", 1)]:
    check(f"stations in tier {t}", float(tiers.get(t, 0)), n, 0)
check("PORI Pescara R", P.loc["Pescara R", "PORI (0-100)"], 100.0, 0.05)
check("PORI Martinsicuro E", P.loc["Martinsicuro E", "PORI (0-100)"], 78.7, 0.1)
check("max AMI gain (PORI points)", ami["PORI_points_avoided"].max(), 9.3, 0.1)
check("robust sentinels", float(sent["Robust_sentinel"].astype(bool).sum()), 2, 0)
print("\nALL CHECKS PASSED" if ok_all else "\nSOME CHECKS FAILED (small differences may come from library versions)")
sys.exit(0 if ok_all else 1)
