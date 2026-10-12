import random
import tkinter as tk
from tkinter import ttk, messagebox

GENETIC_CODE = {}
rows = {
    'UU': ['Phe','Phe','Leu','Leu'], 'UC': ['Ser']*4,
    'UA': ['Tyr','Tyr','Stop','Stop'], 'UG': ['Cys','Cys','Stop','Trp'],
    'CU': ['Leu']*4, 'CC': ['Pro']*4,
    'CA': ['His','His','Gln','Gln'], 'CG': ['Arg']*4,
    'AU': ['Ile','Ile','Ile','Met'], 'AC': ['Thr']*4,
    'AA': ['Asn','Asn','Lys','Lys'], 'AG': ['Ser','Ser','Arg','Arg'],
    'GU': ['Val']*4, 'GC': ['Ala']*4,
    'GA': ['Asp','Asp','Glu','Glu'], 'GG': ['Gly']*4,
}
for prefix, values in rows.items():
    for base, amino in zip('UCAG', values):
        GENETIC_CODE[prefix + base] = amino

class RNAApp(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title('RNA Codon → Amino Acid Translator')
        self.geometry('950x700')
        self.minsize(750, 550)
        self.configure(bg='#f2f5f9')
        style = ttk.Style(self)
        style.theme_use('clam')
        style.configure('TButton', padding=9, font=('Segoe UI', 10))
        style.configure('TLabel', background='#f2f5f9', font=('Segoe UI', 10))
        style.configure('Title.TLabel', font=('Segoe UI', 19, 'bold'))
        style.configure('Treeview', font=('Segoe UI', 10), rowheight=27)
        style.configure('Treeview.Heading', font=('Segoe UI', 10, 'bold'))

        outer = ttk.Frame(self, padding=20)
        outer.pack(fill='both', expand=True)
        ttk.Label(outer, text='RNA → Amino Acid Translator', style='Title.TLabel').pack(anchor='w')
        ttk.Label(outer, text='Generate up to 100 RNA letters and translate consecutive, non-overlapping codons.').pack(anchor='w', pady=(4, 14))

        controls = ttk.Frame(outer)
        controls.pack(fill='x', pady=(0, 12))
        ttk.Label(controls, text='Sequence length:').pack(side='left')
        self.length = tk.IntVar(value=99)
        tk.Spinbox(controls, from_=3, to=100, width=5, textvariable=self.length, font=('Segoe UI', 11)).pack(side='left', padx=(8, 14))
        ttk.Button(controls, text='Generate Random RNA', command=self.generate).pack(side='left', padx=(0, 8))
        ttk.Button(controls, text='Translate', command=self.translate).pack(side='left', padx=(0, 8))
        ttk.Button(controls, text='Copy Amino Acids', command=self.copy_result).pack(side='left')

        ttk.Label(outer, text='RNA sequence (editable; A, U, C, G):').pack(anchor='w')
        self.sequence = tk.Text(outer, height=3, wrap='char', font=('Consolas', 12), relief='solid', borderwidth=1)
        self.sequence.pack(fill='x', pady=(5, 12))

        self.summary = tk.StringVar(value='Generate a sequence to begin.')
        ttk.Label(outer, textvariable=self.summary).pack(anchor='w', pady=(0, 8))
        table_frame = ttk.Frame(outer)
        table_frame.pack(fill='both', expand=True)
        self.table = ttk.Treeview(table_frame, columns=('position','codon','amino'), show='headings')
        for col, heading, width in [('position','Position',110), ('codon','RNA Codon',160), ('amino','Amino Acid',180)]:
            self.table.heading(col, text=heading)
            self.table.column(col, width=width, anchor='center')
        scrollbar = ttk.Scrollbar(table_frame, orient='vertical', command=self.table.yview)
        self.table.configure(yscrollcommand=scrollbar.set)
        self.table.pack(side='left', fill='both', expand=True)
        scrollbar.pack(side='right', fill='y')
        self.table.tag_configure('stop', background='#ffe1df', foreground='#a12121')
        self.table.tag_configure('start', background='#dcfce7', foreground='#166534')

        ttk.Label(outer, text='Amino acid sequence (comma-separated):').pack(anchor='w', pady=(12, 4))
        self.output = tk.Text(outer, height=3, wrap='word', font=('Consolas', 11), relief='solid', borderwidth=1)
        self.output.pack(fill='x')
        self.output.tag_configure('start', foreground='#15803d', font=('Consolas', 11, 'bold'))
        self.output.tag_configure('stop', foreground='#b91c1c', font=('Consolas', 11, 'bold'))
        ttk.Label(outer, text='AUG / Met (START) is green; STOP codons are red. Codons after STOP are also displayed.').pack(anchor='w', pady=(8, 0))
        self.generate()

    def generate(self):
        try:
            length = self.length.get()
            if not 3 <= length <= 100:
                raise ValueError
        except (ValueError, tk.TclError):
            messagebox.showerror('Invalid length', 'Enter a number between 3 and 100.')
            return
        seq = ''.join(random.choice('AUCG') for _ in range(length))
        self.sequence.delete('1.0', 'end')
        self.sequence.insert('1.0', seq)
        self.translate()

    def translate(self):
        seq = ''.join(self.sequence.get('1.0', 'end').upper().split())
        if not seq or any(ch not in 'AUCG' for ch in seq):
            messagebox.showerror('Invalid sequence', 'Use only RNA letters A, U, C, G.')
            return
        if len(seq) > 100:
            messagebox.showerror('Sequence too long', 'The sequence must contain at most 100 letters.')
            return
        for item in self.table.get_children():
            self.table.delete(item)
        acids = []
        for i in range(0, len(seq) - 2, 3):
            codon = seq[i:i+3]
            acid = GENETIC_CODE[codon]
            acids.append(acid)
            self.table.insert('', 'end', values=(f'{i+1}–{i+3}', codon, acid), tags=('stop',) if acid == 'Stop' else ('start',) if acid == 'Met' else ())
        remainder = len(seq) % 3
        self.summary.set(f'{len(seq)} nucleotides  •  {len(acids)} complete codons  •  {sum(a == "Met" for a in acids)} start/Met codons  •  {sum(a == "Stop" for a in acids)} stop codons' + (f'  •  {remainder} trailing letter(s) ignored' if remainder else ''))
        self.output.delete('1.0', 'end')
        for index, acid in enumerate(acids):
            if index:
                self.output.insert('end', ', ')
            label = 'Met (START)' if acid == 'Met' else acid
            tag = 'start' if acid == 'Met' else 'stop' if acid == 'Stop' else None
            self.output.insert('end', label, (tag,) if tag else ())

    def copy_result(self):
        text = self.output.get('1.0', 'end').strip()
        if text:
            self.clipboard_clear()
            self.clipboard_append(text)
            self.update()

if __name__ == '__main__':
    RNAApp().mainloop()
