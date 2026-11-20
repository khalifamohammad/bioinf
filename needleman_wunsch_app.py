"""
Needleman-Wunsch global DNA sequence aligner with a tkinter popup GUI.

Run:   python needleman_wunsch_app.py          -> opens the GUI
       python needleman_wunsch_app.py --cli    -> prints result in the terminal
"""
import sys

S1_DEFAULT = "ACCGTGAAGCCAATAC"
S2_DEFAULT = "AGCGTGCAGCCAATAC"


# --------------------------------------------------------------------------
# Core algorithm
# --------------------------------------------------------------------------
def needleman_wunsch(s1, s2, match=1, mismatch=-1, gap=-2):
    """Return (aligned1, aligned2, score, score_matrix, traceback_path)."""
    n, m = len(s1), len(s2)
    F = [[0] * (m + 1) for _ in range(n + 1)]

    # Initialisation: leading gaps
    for i in range(1, n + 1):
        F[i][0] = i * gap
    for j in range(1, m + 1):
        F[0][j] = j * gap

    # Fill
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            s = match if s1[i - 1] == s2[j - 1] else mismatch
            F[i][j] = max(F[i - 1][j - 1] + s,   # diagonal
                          F[i - 1][j] + gap,     # up   (gap in s2)
                          F[i][j - 1] + gap)     # left (gap in s1)

    # Traceback (prefers diagonal, then up, then left)
    a1, a2, path = [], [], [(n, m)]
    i, j = n, m
    while i > 0 or j > 0:
        if i > 0 and j > 0:
            s = match if s1[i - 1] == s2[j - 1] else mismatch
            if F[i][j] == F[i - 1][j - 1] + s:
                a1.append(s1[i - 1]); a2.append(s2[j - 1])
                i -= 1; j -= 1
                path.append((i, j)); continue
        if i > 0 and F[i][j] == F[i - 1][j] + gap:
            a1.append(s1[i - 1]); a2.append("-")
            i -= 1
        else:
            a1.append("-"); a2.append(s2[j - 1])
            j -= 1
        path.append((i, j))

    return "".join(reversed(a1)), "".join(reversed(a2)), F[n][m], F, path


def format_alignment(a1, a2):
    mid = "".join("|" if x == y else (" " if "-" in (x, y) else "*")
                  for x, y in zip(a1, a2))
    return f"S1: {a1}\n    {mid}\nS2: {a2}"


def stats(a1, a2):
    L = len(a1)
    matches = sum(x == y for x, y in zip(a1, a2))
    gaps = sum("-" in (x, y) for x, y in zip(a1, a2))
    mism = L - matches - gaps
    return L, matches, mism, gaps, 100.0 * matches / L


def clean(seq):
    seq = "".join(seq.split()).upper()
    if not seq or set(seq) - set("ACGT"):
        raise ValueError("Sequences must be non-empty and contain only A, C, G, T.")
    return seq


# --------------------------------------------------------------------------
# GUI
# --------------------------------------------------------------------------
def run_gui():
    import tkinter as tk
    from tkinter import ttk, messagebox

    root = tk.Tk()
    root.title("Needleman-Wunsch DNA Aligner")
    root.geometry("980x720")

    frm = ttk.Frame(root, padding=10)
    frm.pack(fill="both", expand=True)

    ttk.Label(frm, text="Sequence 1:").grid(row=0, column=0, sticky="w")
    e1 = ttk.Entry(frm, width=60); e1.insert(0, S1_DEFAULT)
    e1.grid(row=0, column=1, columnspan=6, sticky="we", pady=2)

    ttk.Label(frm, text="Sequence 2:").grid(row=1, column=0, sticky="w")
    e2 = ttk.Entry(frm, width=60); e2.insert(0, S2_DEFAULT)
    e2.grid(row=1, column=1, columnspan=6, sticky="we", pady=2)

    params = {}
    for k, (name, val) in enumerate([("Match", 1), ("Mismatch", -1), ("Gap", -2)]):
        ttk.Label(frm, text=name + ":").grid(row=2, column=2 * k, sticky="e", pady=6)
        v = tk.StringVar(value=str(val))
        ttk.Spinbox(frm, from_=-20, to=20, width=5, textvariable=v)\
            .grid(row=2, column=2 * k + 1, sticky="w")
        params[name] = v

    out = tk.Text(frm, height=8, font=("Courier New", 12))
    out.grid(row=4, column=0, columnspan=7, sticky="we", pady=6)

    ttk.Label(frm, text="Score matrix (red = traceback path):")\
        .grid(row=5, column=0, columnspan=7, sticky="w")
    canvas_frame = ttk.Frame(frm)
    canvas_frame.grid(row=6, column=0, columnspan=7, sticky="nsew")
    canvas = tk.Canvas(canvas_frame, bg="white")
    hs = ttk.Scrollbar(canvas_frame, orient="horizontal", command=canvas.xview)
    vs = ttk.Scrollbar(canvas_frame, orient="vertical", command=canvas.yview)
    canvas.configure(xscrollcommand=hs.set, yscrollcommand=vs.set)
    hs.pack(side="bottom", fill="x"); vs.pack(side="right", fill="y")
    canvas.pack(side="left", fill="both", expand=True)

    frm.columnconfigure(6, weight=1)
    frm.rowconfigure(6, weight=1)

    def draw_matrix(s1, s2, F, path):
        canvas.delete("all")
        cell = 34
        on_path = set(path)
        heads2 = " " + s2
        for j, ch in enumerate(heads2):
            canvas.create_text((j + 1.5) * cell, cell / 2, text=ch, font=("Arial", 11, "bold"))
        for i in range(len(s1) + 1):
            lab = " " if i == 0 else s1[i - 1]
            canvas.create_text(cell / 2, (i + 1.5) * cell, text=lab, font=("Arial", 11, "bold"))
            for j in range(len(s2) + 1):
                x, y = (j + 1) * cell, (i + 1) * cell
                hit = (i, j) in on_path
                canvas.create_rectangle(x, y, x + cell, y + cell,
                                        fill="#ffd0d0" if hit else "white", outline="#999")
                canvas.create_text(x + cell / 2, y + cell / 2, text=str(F[i][j]),
                                   fill="red" if hit else "black", font=("Arial", 10))
        canvas.configure(scrollregion=canvas.bbox("all"))

    def align():
        try:
            s1, s2 = clean(e1.get()), clean(e2.get())
            ma, mi, ga = (int(params[k].get()) for k in ("Match", "Mismatch", "Gap"))
        except ValueError as err:
            messagebox.showerror("Invalid input", str(err))
            return
        a1, a2, score, F, path = needleman_wunsch(s1, s2, ma, mi, ga)
        L, mt, mm, gp, ident = stats(a1, a2)
        report = (f"{format_alignment(a1, a2)}\n\n"
                  f"Score: {score}   Length: {L}   Matches: {mt}   "
                  f"Mismatches: {mm}   Gaps: {gp}   Identity: {ident:.1f}%")
        out.delete("1.0", "end"); out.insert("end", report)
        draw_matrix(s1, s2, F, path)

        # Popup window with the final result
        messagebox.showinfo("Alignment result", report)

    ttk.Button(frm, text="Align", command=align)\
        .grid(row=3, column=0, columnspan=7, pady=6)

    root.mainloop()


def run_cli():
    a1, a2, score, *_ = needleman_wunsch(S1_DEFAULT, S2_DEFAULT)
    L, mt, mm, gp, ident = stats(a1, a2)
    print(format_alignment(a1, a2))
    print(f"\nScore={score} Length={L} Matches={mt} Mismatches={mm} Gaps={gp} Identity={ident:.1f}%")


if __name__ == "__main__":
    if "--cli" in sys.argv:
        run_cli()
    else:
        try:
            run_gui()
        except Exception as exc:  # no display / tkinter missing
            print(f"GUI unavailable ({exc}); falling back to CLI.\n")
            run_cli()
