"""
DNA pattern generator for gene promoters (PromKappa-style).

For every sliding window of a sequence we compute
    x = C+G %                      (cg_percent)
    y = Kappa Index of Coincidence (kappa_ic)  <- the "information content"
The cloud of (x, y) points is the DNA pattern; its mean is the "center of weight".

Requires: matplotlib  (pip install matplotlib).  GUI uses tkinter (ships with Python).

  python dna_pattern_app.py            -> GUI
  python dna_pattern_app.py --cli      -> self-test + saves PNG charts (no GUI needed)
  python dna_pattern_app.py --cli promoters.fa   -> also builds patterns from a FASTA file
"""
import sys
import urllib.request

# Test sequence of the assignment: CGGACTGATCTATCT + (A x 27) + CGTAGCATCTATCGATCTATCTCTAGCGATCTA
TEST_SEQ = "CGGACTGATCTATCT" + "A" * 27 + "CGTAGCATCTATCGATCTATCTCTAGCGATCTA"
TEST_WINDOW = 30
FASTA_URL = ("https://raw.githubusercontent.com/Gagniuc/"
             "PromKappa-2.0-makes-Objective-Digital-Stains/main/gene_promoters_complete.fa")


# --------------------------------------------------------------------------
# Core functions
# --------------------------------------------------------------------------
def cg_percent(seq):
    """C+G percentage of a sequence."""
    seq = seq.upper()
    return 100.0 * (seq.count("C") + seq.count("G")) / len(seq)


def kappa_ic(seq):
    """Kappa Index of Coincidence, a port of PromKappa's IC() function (Kappa.bas).

    For every shift u = 1..L-1 the sequence is compared with itself shifted by u;
    the percentage of coincident positions is accumulated and averaged over L-1 shifts.
    """
    seq = seq.upper()
    L = len(seq)
    if L < 2:
        return 0.0
    total = 0.0
    for u in range(1, L):
        coincidences = sum(seq[i] == seq[i + u] for i in range(L - u))
        total += coincidences / (L - u) * 100.0
    return round(total / (L - 1), 2)


def sliding_pattern(seq, window=30, step=1):
    """Return the pattern: list of (C+G%, Kappa IC), one per window.
    Like PromKappa, windows start at 1 .. L-window (the very last window is not used)."""
    seq = "".join(seq.split()).upper()
    if window >= len(seq):
        return [(cg_percent(seq), kappa_ic(seq))]
    return [(cg_percent(seq[i:i + window]), kappa_ic(seq[i:i + window]))
            for i in range(0, len(seq) - window, step)]


def center_of_weight(pattern):
    """Centroid (mean x, mean y) of the pattern points."""
    n = len(pattern)
    return (sum(p[0] for p in pattern) / n, sum(p[1] for p in pattern) / n)


def read_fasta(text):
    records, name, buf = [], None, []
    for line in text.splitlines():
        line = line.strip()
        if line.startswith(">"):
            if name is not None:
                records.append((name, "".join(buf).upper()))
            name, buf = line[1:], []
        elif line:
            buf.append(line)
    if name is not None:
        records.append((name, "".join(buf).upper()))
    return records


def load_fasta(source):
    """source = local path or http(s) URL."""
    if source.startswith("http"):
        with urllib.request.urlopen(source, timeout=60) as r:
            return read_fasta(r.read().decode("latin-1"))
    with open(source, encoding="latin-1") as f:
        return read_fasta(f.read())


# --------------------------------------------------------------------------
# Charts (matplotlib figures, used by both GUI and CLI)
# --------------------------------------------------------------------------
def plot_pattern(ax, pattern, title="DNA pattern"):
    xs, ys = zip(*pattern)
    cx, cy = center_of_weight(pattern)
    ax.scatter(xs, ys, s=18, alpha=0.7, label="window")
    ax.plot(xs, ys, lw=0.5, color="gray", alpha=0.5)
    ax.scatter([cx], [cy], s=120, c="black", marker="o", label=f"center ({cx:.2f}, {cy:.2f})")
    ax.set_xlabel("C+G %"); ax.set_ylabel("Kappa IC (information content)")
    ax.set_title(title); ax.grid(alpha=0.3); ax.legend(fontsize=8)


def plot_centers(ax, centers, labels=None, title="Centers of weight"):
    xs, ys = zip(*centers)
    ax.scatter(xs, ys, s=25, c="crimson")
    if labels and len(centers) <= 15:
        for (x, y), t in zip(centers, labels):
            ax.annotate(t, (x, y), fontsize=7, xytext=(3, 3), textcoords="offset points")
    ax.set_xlabel("Average C+G %"); ax.set_ylabel("Average Kappa IC")
    ax.set_title(title); ax.grid(alpha=0.3)


# --------------------------------------------------------------------------
# Self-test on the sequence given in the assignment
# --------------------------------------------------------------------------
def self_test():
    S = TEST_SEQ
    first = S[:TEST_WINDOW]
    pat = sliding_pattern(S, TEST_WINDOW)
    lines = [f"Sequence length: {len(S)}  | window: {TEST_WINDOW}  | windows in pattern: {len(pat)}",
             f"First window  : C+G% = {cg_percent(first):.2f}   Kappa IC = {kappa_ic(first):.2f}",
             f"Whole sequence: C+G% = {cg_percent(S):.2f}   Kappa IC = {kappa_ic(S):.2f}",
             "Reference values given in the assignment: CG ~ 29.27, Ic ~ 27.53",
             "Center of weight of the pattern: (%.2f, %.2f)" % center_of_weight(pat)]
    return "\n".join(lines)


# --------------------------------------------------------------------------
# GUI
# --------------------------------------------------------------------------
def run_gui():
    import tkinter as tk
    from tkinter import ttk, filedialog, messagebox
    from matplotlib.figure import Figure
    from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

    root = tk.Tk()
    root.title("DNA Pattern Generator (C+G% vs Kappa IC)")
    root.geometry("1100x800")
    state = {"records": []}

    top = ttk.Frame(root, padding=8); top.pack(fill="x")
    ttk.Label(top, text="Sequence:").grid(row=0, column=0, sticky="w")
    seq_txt = tk.Text(top, height=3, width=90); seq_txt.insert("1.0", TEST_SEQ)
    seq_txt.grid(row=0, column=1, columnspan=6, sticky="we")
    ttk.Label(top, text="Window:").grid(row=1, column=0, sticky="e")
    win_v = tk.StringVar(value=str(TEST_WINDOW))
    ttk.Spinbox(top, from_=5, to=1000, width=6, textvariable=win_v).grid(row=1, column=1, sticky="w")
    ttk.Label(top, text="Step:").grid(row=1, column=2, sticky="e")
    step_v = tk.StringVar(value="1")
    ttk.Spinbox(top, from_=1, to=100, width=6, textvariable=step_v).grid(row=1, column=3, sticky="w")

    nb = ttk.Notebook(root); nb.pack(fill="both", expand=True)
    fig1, fig2 = Figure(figsize=(6, 4.5)), Figure(figsize=(6, 4.5))
    tabs = []
    for fig, name in ((fig1, "Pattern"), (fig2, "Centers of weight")):
        f = ttk.Frame(nb); nb.add(f, text=name)
        c = FigureCanvasTkAgg(fig, master=f); c.get_tk_widget().pack(fill="both", expand=True)
        tabs.append(c)
    log = tk.Text(root, height=7, font=("Courier New", 10)); log.pack(fill="x")

    def params():
        return int(win_v.get()), int(step_v.get())

    def make_pattern():
        try:
            seq = "".join(seq_txt.get("1.0", "end").split()).upper()
            if not seq or set(seq) - set("ACGT"):
                raise ValueError("Sequence must contain only A, C, G, T.")
            w, st = params()
            pat = sliding_pattern(seq, w, st)
        except ValueError as e:
            messagebox.showerror("Invalid input", str(e)); return
        fig1.clear(); plot_pattern(fig1.add_subplot(111), pat, f"Pattern (window={w})")
        tabs[0].draw(); nb.select(0)
        cx, cy = center_of_weight(pat)
        log.delete("1.0", "end")
        log.insert("end", f"Windows: {len(pat)}\nCenter of weight: C+G%={cx:.2f}, Kappa IC={cy:.2f}\n")
        # Add this center to the second chart
        state.setdefault("centers", []).append((cx, cy))
        fig2.clear(); plot_centers(fig2.add_subplot(111), state["centers"],
                                   [str(i + 1) for i in range(len(state["centers"]))])
        tabs[1].draw()

    def run_selftest():
        text = self_test()
        log.delete("1.0", "end"); log.insert("end", text)
        messagebox.showinfo("Self-test", text)

    def load_promoters(url=False):
        try:
            src = FASTA_URL if url else filedialog.askopenfilename(
                title="Open promoter FASTA (e.g. gene_promoters_complete.fa)")
            if not src:
                return
            state["records"] = load_fasta(src)
        except Exception as e:
            messagebox.showerror("Load error", str(e)); return
        log.delete("1.0", "end")
        log.insert("end", f"Loaded {len(state['records'])} promoters. First: {state['records'][0][0]}\n")

    def use_first_promoter():
        if not state["records"]:
            messagebox.showinfo("No data", "Load a promoter FASTA first."); return
        name, s = state["records"][0]
        seq_txt.delete("1.0", "end"); seq_txt.insert("1.0", s)
        win_v.set("100")
        make_pattern()
        log.insert("end", f"Promoter: {name}\n")

    def centers_of_many():
        if not state["records"]:
            messagebox.showinfo("No data", "Load a promoter FASTA first."); return
        w, st = params()
        recs = state["records"][:200]
        centers, labels = [], []
        for name, s in recs:
            if len(s) > w:
                centers.append(center_of_weight(sliding_pattern(s, w, st)))
                labels.append(name.split(";")[0].split(")")[-1].strip())
        state["centers"] = centers
        fig2.clear(); plot_centers(fig2.add_subplot(111), centers, labels,
                                   f"Centers of weight of {len(centers)} promoters")
        tabs[1].draw(); nb.select(1)

    bar = ttk.Frame(top); bar.grid(row=2, column=0, columnspan=7, pady=6, sticky="w")
    for text, cmd in (("Generate pattern", make_pattern), ("Self-test (assignment values)", run_selftest),
                      ("Load FASTA file…", lambda: load_promoters(False)),
                      ("Download PromKappa promoters", lambda: load_promoters(True)),
                      ("Pattern of first promoter", use_first_promoter),
                      ("Centers of 200 promoters", centers_of_many)):
        ttk.Button(bar, text=text, command=cmd).pack(side="left", padx=3)

    make_pattern()
    root.mainloop()


# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------
def run_cli(fasta=None):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    print(self_test())
    pat = sliding_pattern(TEST_SEQ, TEST_WINDOW)
    fig, ax = plt.subplots(figsize=(6, 4.5))
    plot_pattern(ax, pat, "Pattern of test sequence (window=30)")
    fig.savefig("pattern_test_sequence.png", dpi=120)

    if fasta:
        recs = load_fasta(fasta)
        print(f"\nLoaded {len(recs)} promoters from {fasta}")
        name, s = recs[0]
        p = sliding_pattern(s, 100)
        fig, ax = plt.subplots(figsize=(6, 4.5))
        plot_pattern(ax, p, f"Promoter pattern: {name[:40]} (window=100)")
        fig.savefig("pattern_promoter.png", dpi=120)
        centers = [center_of_weight(sliding_pattern(x, 100)) for _, x in recs[:200] if len(x) > 100]
        fig, ax = plt.subplots(figsize=(6, 4.5))
        plot_centers(ax, centers, title=f"Centers of weight of {len(centers)} promoters")
        fig.savefig("centers.png", dpi=120)
        print("Center of weight of first promoter: (%.2f, %.2f)" % center_of_weight(p))
    print("Charts saved in the current folder.")


if __name__ == "__main__":
    if "--cli" in sys.argv:
        args = [a for a in sys.argv[1:] if a != "--cli"]
        run_cli(args[0] if args else None)
    else:
        try:
            run_gui()
        except Exception as exc:
            print(f"GUI unavailable ({exc}); running CLI instead.\n")
            run_cli()
