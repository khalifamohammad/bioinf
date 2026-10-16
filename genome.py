"""Interactive whole-genome overlapping codon analyzer. Python 3 standard library only."""
import gzip
import tkinter as tk
from tkinter import ttk, filedialog, messagebox
from collections import Counter
from pathlib import Path

CODE_ROWS = [
('TTT TTC TTA TTG','Phe Phe Leu Leu'),('TCT TCC TCA TCG','Ser Ser Ser Ser'),('TAT TAC TAA TAG','Tyr Tyr Stop Stop'),('TGT TGC TGA TGG','Cys Cys Stop Trp'),
('CTT CTC CTA CTG','Leu Leu Leu Leu'),('CCT CCC CCA CCG','Pro Pro Pro Pro'),('CAT CAC CAA CAG','His His Gln Gln'),('CGT CGC CGA CGG','Arg Arg Arg Arg'),
('ATT ATC ATA ATG','Ile Ile Ile Met'),('ACT ACC ACA ACG','Thr Thr Thr Thr'),('AAT AAC AAA AAG','Asn Asn Lys Lys'),('AGT AGC AGA AGG','Ser Ser Arg Arg'),
('GTT GTC GTA GTG','Val Val Val Val'),('GCT GCC GCA GCG','Ala Ala Ala Ala'),('GAT GAC GAA GAG','Asp Asp Glu Glu'),('GGT GGC GGA GGG','Gly Gly Gly Gly')]
GENETIC_CODE = {c:a for cs,aas in CODE_ROWS for c,a in zip(cs.split(),aas.split())}
FOODS = [
    'Cooking oils: essentially no protein or amino acids.',
    'Table sugar: essentially no protein or amino acids.',
    'Most fresh fruits: generally low total protein per serving.',
    'Refined starches: often relatively low in protein, depending on the food.'
]

def read_genome(path):
    opener = gzip.open if str(path).lower().endswith('.gz') else open
    records, current = [], []
    with opener(path, 'rt', encoding='utf-8-sig') as handle:
        for line in handle:
            line = line.strip()
            if not line: continue
            if line.startswith('>'):
                if current: records.append(''.join(current)); current = []
                continue
            if line.startswith(';'): continue
            current.append(''.join(ch for ch in line.upper() if ch.isalpha()))
    if current: records.append(''.join(current))
    return 'NNN'.join(records)

class GenomeApp:
    def __init__(self, root):
        self.root = root
        root.title('Genome Analyzer | Live Sliding Window')
        root.geometry('1040x790')
        root.minsize(790,620)
        self.sequence = ''
        self.path = ''
        self.counts = Counter()
        self.skipped = 0
        self.position = 0
        self.running = False
        self.after_id = None
        self.speed = tk.IntVar(value=10)
        self.path_var = tk.StringVar(value='No genome loaded')
        self.info_var = tk.StringVar(value='Load a FASTA, FNA, TXT or .gz genome to begin.')
        self.codon_var = tk.StringVar(value='Codon: —     Amino acid: —')
        self.position_var = tk.StringVar(value='Window: 0 / 0')
        self.build()
        root.protocol('WM_DELETE_WINDOW', self.close)

    def build(self):
        top = ttk.Frame(self.root, padding=12)
        top.pack(fill='x')
        ttk.Label(top,text='LIVE GENOME CODON ANALYZER',font=('Segoe UI',18,'bold')).pack(anchor='w')
        bar = ttk.Frame(top)
        bar.pack(fill='x',pady=(8,0))
        ttk.Button(bar,text='Browse genome…',command=self.browse).pack(side='left')
        ttk.Label(bar,textvariable=self.path_var).pack(side='left',padx=12)
        ttk.Label(top,textvariable=self.info_var,foreground='#52677c').pack(anchor='w',pady=(8,0))

        viewer = ttk.LabelFrame(self.root,text='Sliding window — step size 1 nucleotide',padding=12)
        viewer.pack(fill='x',padx=12,pady=(0,10))
        ttk.Label(viewer,text='The highlighted 3 bases are the current overlapping codon.').pack(anchor='w')
        self.sequence_text = tk.Text(viewer,height=3,wrap='char',font=('Consolas',16),state='disabled',bg='#111d2b',fg='#d9e5f2',relief='flat',padx=12,pady=9)
        self.sequence_text.pack(fill='x',pady=7)
        self.sequence_text.tag_configure('current',background='#16a085',foreground='white',font=('Consolas',16,'bold'))
        self.sequence_text.tag_configure('invalid',background='#b03a2e',foreground='white')
        ttk.Label(viewer,textvariable=self.codon_var,font=('Segoe UI',12,'bold')).pack(anchor='w')
        ttk.Label(viewer,textvariable=self.position_var).pack(anchor='w',pady=(3,0))
        self.progress = ttk.Progressbar(viewer,mode='determinate',maximum=100)
        self.progress.pack(fill='x',pady=8)
        controls = ttk.Frame(viewer)
        controls.pack(fill='x')
        self.play_button = ttk.Button(controls,text='▶ Play',command=self.play)
        self.play_button.pack(side='left',padx=(0,5))
        ttk.Button(controls,text='⏸ Pause',command=self.pause).pack(side='left',padx=5)
        ttk.Button(controls,text='↺ Restart',command=self.restart).pack(side='left',padx=5)
        ttk.Button(controls,text='⏭ Analyze all',command=self.analyze_all).pack(side='left',padx=5)
        ttk.Label(controls,text='Windows per frame:').pack(side='left',padx=(20,6))
        ttk.Combobox(controls,textvariable=self.speed,values=(1,10,100,1000,10000),width=7,state='readonly').pack(side='left')
        ttk.Label(viewer,text='Jump to window (1-based):').pack(anchor='w',pady=(12,0))
        jump = ttk.Frame(viewer)
        jump.pack(fill='x',pady=(4,0))
        self.seek = ttk.Scale(jump,from_=1,to=2,orient='horizontal',command=self.preview_seek)
        self.seek.pack(side='left',fill='x',expand=True)
        ttk.Button(jump,text='Jump & scan to here',command=self.jump).pack(side='left',padx=(12,0))
        ttk.Label(viewer,text='Jump rescans from the beginning to preserve accurate cumulative frequencies.',foreground='#52677c').pack(anchor='w',pady=(4,0))

        body = ttk.Frame(self.root,padding=(12,0,12,12))
        body.pack(fill='both',expand=True)
        left = ttk.LabelFrame(body,text='Top 10 codons (cumulative counts)',padding=8)
        left.pack(side='left',fill='both',expand=True,padx=(0,6))
        columns = ('rank','codon','count','freq','amino')
        self.table = ttk.Treeview(left,columns=columns,show='headings',height=12)
        for key,label,width in [('rank','#',40),('codon','Codon',85),('count','Count',95),('freq','Frequency',95),('amino','Amino acid',100)]:
            self.table.heading(key,text=label)
            self.table.column(key,width=width,anchor='center',stretch=True)
        self.table.tag_configure('start',foreground='#16834a')
        self.table.tag_configure('stop',foreground='#c0392b')
        self.table.pack(fill='both',expand=True)
        right = ttk.LabelFrame(body,text='Amino acids & illustrative food comparisons',padding=8)
        right.pack(side='left',fill='both',expand=True,padx=(6,0))
        self.notes = tk.Text(right,wrap='word',font=('Segoe UI',10),height=14)
        self.notes.pack(fill='both',expand=True)
        ttk.Button(right,text='Save results (.txt)',command=self.save_results).pack(anchor='e',pady=(8,0))
        self.refresh_notes()

    def browse(self):
        path = filedialog.askopenfilename(filetypes=[('Genome files','*.fasta *.fa *.fna *.txt *.gz'),('All files','*.*')])
        if not path: return
        self.pause()
        try: seq = read_genome(path)
        except (OSError,UnicodeError) as exc:
            messagebox.showerror('Could not load genome',str(exc)); return
        self.sequence,self.path = seq,path
        self.path_var.set(Path(path).name)
        self.seek.configure(to=max(1,len(seq)-2))
        self.restart()
        self.info_var.set(f'{len(seq):,} bases loaded | {max(0,len(seq)-2):,} overlapping windows')
        if len(seq)<3: messagebox.showwarning('Short sequence','Genome must contain at least 3 bases.')

    def pause(self):
        self.running = False
        if self.after_id is not None:
            self.root.after_cancel(self.after_id)
            self.after_id = None

    def restart(self):
        self.pause()
        self.position = 0
        self.counts.clear()
        self.skipped = 0
        self.update_display()

    def play(self):
        if len(self.sequence)<3:
            messagebox.showinfo('Load genome','Load a genome with at least 3 bases first.')
            return
        if self.position>=len(self.sequence)-2: self.restart()
        if not self.running:
            self.running = True
            self.tick()

    def scan_until(self,target):
        limit = max(0,len(self.sequence)-2)
        target = min(target,limit)
        while self.position < target:
            codon = self.sequence[self.position:self.position+3]
            if len(codon)==3 and all(c in 'ACGT' for c in codon): self.counts[codon]+=1
            else: self.skipped+=1
            self.position+=1

    def tick(self):
        if not self.running: return
        self.scan_until(self.position+int(self.speed.get()))
        self.update_display()
        if self.position>=max(0,len(self.sequence)-2):
            self.running=False
            self.after_id=None
        else: self.after_id=self.root.after(45,self.tick)

    def analyze_all(self):
        if len(self.sequence)<3: return
        self.pause()
        self.scan_until(len(self.sequence)-2)
        self.update_display()

    def preview_seek(self,value):
        if not self.sequence or self.running: return
        i = max(0,min(int(float(value))-1,len(self.sequence)-3))
        self.show_window(i,preview=True)

    def jump(self):
        if len(self.sequence)<3: return
        self.restart()
        target = int(self.seek.get())
        self.scan_until(target)
        self.update_display()

    def show_window(self,i,preview=False):
        if len(self.sequence)<3:
            self.render_sequence('')
            self.codon_var.set('Codon: —     Amino acid: —')
            return
        i = min(max(0,i),len(self.sequence)-3)
        codon = self.sequence[i:i+3]
        valid = all(c in 'ACGT' for c in codon)
        amino = GENETIC_CODE[codon] if valid else 'Invalid / skipped'
        start = max(0,i-45)
        end = min(len(self.sequence),i+48)
        fragment = self.sequence[start:end]
        self.render_sequence(fragment,i-start,valid)
        self.codon_var.set(f'{"Preview — " if preview else ""}Codon: {codon}    →    {amino}' + (' (START / Met)' if codon=='ATG' else ' (STOP)' if amino=='Stop' else ''))

    def render_sequence(self,fragment,offset=None,valid=True):
        self.sequence_text.configure(state='normal')
        self.sequence_text.delete('1.0','end')
        self.sequence_text.insert('1.0',fragment)
        if offset is not None:
            self.sequence_text.tag_add('current' if valid else 'invalid',f'1.0+{offset}c',f'1.0+{offset+3}c')
        self.sequence_text.configure(state='disabled')

    def update_display(self):
        total_windows = max(0,len(self.sequence)-2)
        self.show_window(max(0,self.position-1))
        self.position_var.set(f'Windows analyzed: {self.position:,} / {total_windows:,}   |   Valid: {sum(self.counts.values()):,}   |   Skipped: {self.skipped:,}')
        self.progress['value'] = self.position/total_windows*100 if total_windows else 0
        if not self.running and total_windows: self.seek.set(max(1,self.position))
        for item in self.table.get_children(): self.table.delete(item)
        valid = sum(self.counts.values())
        top = sorted(self.counts.items(),key=lambda item:(-item[1],item[0]))[:10]
        for rank,(codon,count) in enumerate(top,1):
            aa = GENETIC_CODE[codon]
            self.table.insert('','end',values=(rank,codon,f'{count:,}',f'{100*count/valid:.2f}%',aa),tags=('start',) if codon=='ATG' else ('stop',) if aa=='Stop' else ())
        self.refresh_notes(top)

    def refresh_notes(self,top=()):
        acids = list(dict.fromkeys(GENETIC_CODE[c] for c,_ in top if GENETIC_CODE[c]!='Stop'))
        question = ('Given these amino acids that appear most frequently according to the codon counts: '
                    + ', '.join(acids) + ', which common foods are relatively low in or lack these amino acids?')
        content = 'AMINO ACIDS FROM TOP CODONS\n' + (', '.join(acids) if acids else 'None counted yet')
        content += '\n\nSIMULATED AI QUESTION\n' + question
        content += '\n\nILLUSTRATIVE FOOD EXAMPLES\n' + '\n'.join('• '+s for s in FOODS)
        content += ('\n\nIMPORTANT: Overlapping triplets across a genome do not necessarily correspond '
                    'to translated protein codons. Their frequencies do not measure dietary amino-acid needs. '
                    'The food examples are general, not amino-acid-specific nutritional measurements.')
        self.notes.configure(state='normal')
        self.notes.delete('1.0','end')
        self.notes.insert('1.0',content)
        self.notes.configure(state='disabled')

    def save_results(self):
        if not self.sequence:
            messagebox.showinfo('No genome','Load and analyze a genome first.'); return
        path = filedialog.asksaveasfilename(defaultextension='.txt',filetypes=[('Text report','*.txt')])
        if not path: return
        total=sum(self.counts.values())
        top=sorted(self.counts.items(),key=lambda item:(-item[1],item[0]))[:10]
        lines=[f'Genome: {self.path}',f'Bases: {len(self.sequence):,}',f'Windows analyzed: {self.position:,}',f'Valid: {total:,}',f'Skipped: {self.skipped:,}','', 'Top 10 codons:']
        lines += [f'{i:2}. {c}  {n:,}  {100*n/total:.2f}%  {GENETIC_CODE[c]}' for i,(c,n) in enumerate(top,1)] if total else []
        lines += ['','',self.notes.get('1.0','end').strip()]
        try: Path(path).write_text('\n'.join(lines),encoding='utf-8')
        except OSError as exc: messagebox.showerror('Save error',str(exc)); return
        messagebox.showinfo('Saved',f'Report saved to:\n{path}')

    def close(self):
        self.pause()
        self.root.destroy()

if __name__=='__main__':
    root=tk.Tk()
    GenomeApp(root)
    root.mainloop()
