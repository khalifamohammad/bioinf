"""DNA Reconstruction Lab — standard-library Tkinter desktop application."""
import gzip
import random
import threading
import queue
import time
from pathlib import Path
import tkinter as tk
from tkinter import ttk, filedialog, messagebox

BASE = Path(__file__).resolve().parent
DEFAULT = next((BASE / 'data').glob('*.gz'), None)
COLORS = {'A':'#46d4a6','C':'#64a9ff','G':'#f5c56a','T':'#f18aab','N':'#e6e6e6'}

def read_fasta(path):
    opener = gzip.open if str(path).lower().endswith('.gz') else open
    with opener(path, 'rt', encoding='utf-8-sig') as f:
        records=[]; current=[]; name=''
        for line in f:
            line=line.strip()
            if line.startswith('>'):
                if current: records.append((name, ''.join(current))); current=[]
                name=line[1:]
            elif line and not line.startswith(';'):
                current.append(line.upper())
        if current: records.append((name, ''.join(current)))
    if not records: raise ValueError('No FASTA sequence found')
    name, seq=records[0]
    if set(seq)-set('ACGTN'): raise ValueError('Genome contains unsupported symbols')
    return name, seq, len(records)

def make_reads(seq, count, low, high, seed):
    rng=random.Random(seed)
    reads=[]
    for i in range(count):
        length=rng.randint(low, min(high,len(seq)))
        start=rng.randint(0,len(seq)-length)
        reads.append({'id':i+1,'sequence':seq[start:start+length],'start':start})
    rng.shuffle(reads)
    return reads

def overlap(a,b,minimum):
    for n in range(min(len(a),len(b)),minimum-1,-1):
        if a.endswith(b[:n]): return n
    return 0

def assemble(reads, minimum, emit):
    """Greedy exact-overlap assembly. Reads' source coordinates are never used."""
    pieces=list(dict.fromkeys(r['sequence'] for r in reads))
    # Remove reads wholly contained in longer reads, independent of original position.
    pieces.sort(key=len,reverse=True)
    kept=[]
    for p in pieces:
        if not any(p in longer for longer in kept): kept.append(p)
    pieces=kept
    merges=0; start=time.perf_counter()
    while len(pieces)>1:
        best=(0,-1,-1)
        for i,a in enumerate(pieces):
            for j,b in enumerate(pieces):
                if i==j: continue
                score=overlap(a,b,minimum)
                if score>best[0]: best=(score,i,j)
        n,i,j=best
        if i<0: break
        a,b=pieces[i],pieces[j]
        combined=a+b[n:]
        pieces=[p for k,p in enumerate(pieces) if k!=i and k!=j]
        pieces.append(combined)
        # Eliminate newly contained contigs.
        pieces=[p for p in pieces if p==combined or p not in combined]
        merges+=1
        emit(('merge', merges, a, b, n, combined, len(pieces), time.perf_counter()-start))
        if merges%10==0: emit(('progress',merges,len(pieces),time.perf_counter()-start))
    return sorted(pieces,key=len,reverse=True),merges,time.perf_counter()-start

class App:
    def __init__(self,root):
        self.root=root; root.title('GenomeLab | DNA Reconstruction'); root.geometry('1100x820'); root.minsize(850,650)
        self.name='';self.genome='';self.seq='';self.reads=[];self.contigs=[];self.busy=False;self.messages=queue.Queue();self.history=[];self.play_index=0;self.playing=False;self.play_after=None
        self.start_var=tk.IntVar(value=1);self.window_var=tk.IntVar(value=120)
        self.count_var=tk.IntVar(value=2000);self.low_var=tk.IntVar(value=100);self.high_var=tk.IntVar(value=120)
        self.seed_var=tk.IntVar(value=42);self.overlap_var=tk.IntVar(value=35)
        self.info=tk.StringVar(value='Load a genome to begin');self.status=tk.StringVar(value='Ready')
        self.build();root.after(120,self.poll)
        if DEFAULT: self.load(DEFAULT)

    def build(self):
        style=ttk.Style();style.theme_use('clam')
        style.configure('TButton',padding=7);style.configure('TNotebook.Tab',padding=(18,9))
        head=ttk.Frame(self.root,padding=15);head.pack(fill='x')
        ttk.Label(head,text='GENOMELAB',font=('Segoe UI',23,'bold')).pack(side='left')
        ttk.Label(head,text='  DNA Reconstruction Studio',font=('Segoe UI',12)).pack(side='left',padx=12)
        ttk.Button(head,text='Open genome…',command=self.browse).pack(side='right')
        ttk.Label(self.root,textvariable=self.info,foreground='#3c7194',padding=(16,0)).pack(anchor='w')
        self.tabs=ttk.Notebook(self.root);self.tabs.pack(fill='both',expand=True,padx=12,pady=12)
        self.genome_tab=ttk.Frame(self.tabs,padding=14);self.samples_tab=ttk.Frame(self.tabs,padding=14);self.assembly_tab=ttk.Frame(self.tabs,padding=14)
        self.tabs.add(self.genome_tab,text='01  Genome Viewer');self.tabs.add(self.samples_tab,text='02  Random Samples');self.tabs.add(self.assembly_tab,text='03  Reconstruction')
        self.live_tab=ttk.Frame(self.tabs,padding=14);self.tabs.add(self.live_tab,text='04  Live Reconstruction')
        self.build_viewer();self.build_samples();self.build_assembly();self.build_live()
        footer=ttk.Frame(self.root,padding=(15,0,15,12));footer.pack(fill='x')
        ttk.Label(footer,textvariable=self.status).pack(side='left')
        self.progress=ttk.Progressbar(footer,mode='indeterminate',length=170);self.progress.pack(side='right')

    def build_viewer(self):
        p=self.genome_tab
        ttk.Label(p,text='Original DNA sequence',font=('Segoe UI',16,'bold')).pack(anchor='w')
        ttk.Label(p,text='Move the slider to inspect the genome. The selected window is highlighted in color.').pack(anchor='w',pady=(4,12))
        controls=ttk.Frame(p);controls.pack(fill='x')
        ttk.Label(controls,text='Window length').pack(side='left')
        ttk.Spinbox(controls,from_=20,to=400,textvariable=self.window_var,width=7,command=self.draw_sequence).pack(side='left',padx=8)
        ttk.Button(controls,text='Use first 2,000 bases',command=lambda:self.select_region(False)).pack(side='left',padx=7)
        ttk.Button(controls,text='Use entire genome',command=lambda:self.select_region(True)).pack(side='left')
        self.region_label=ttk.Label(p,text='No sequence selected');self.region_label.pack(anchor='w',pady=(12,5))
        self.slider=ttk.Scale(p,from_=1,to=2,variable=self.start_var,command=lambda _:self.draw_sequence());self.slider.pack(fill='x',pady=(0,12))
        self.dna_text=tk.Text(p,height=10,wrap='char',font=('Consolas',15),bg='#111d2b',fg='white',padx=16,pady=16)
        self.dna_text.pack(fill='both',expand=True)
        for b,c in COLORS.items():self.dna_text.tag_configure(b,foreground=c)
        self.dna_text.configure(state='disabled')
        ttk.Label(p,text='A = green    C = blue    G = gold    T = pink',foreground='#5a6a7a').pack(anchor='w',pady=10)

    def build_samples(self):
        p=self.samples_tab
        ttk.Label(p,text='Generate random fragments',font=('Segoe UI',16,'bold')).pack(anchor='w')
        ttk.Label(p,text='Samples are randomly drawn from the selected sequence. Start positions are hidden from the assembler.').pack(anchor='w',pady=(4,12))
        bar=ttk.Frame(p);bar.pack(fill='x')
        for label,var in [('Samples',self.count_var),('Min bases',self.low_var),('Max bases',self.high_var),('Random seed',self.seed_var)]:
            item=ttk.Frame(bar);item.pack(side='left',padx=(0,15));ttk.Label(item,text=label).pack(anchor='w');ttk.Entry(item,textvariable=var,width=10).pack()
        ttk.Button(bar,text='Generate samples',command=self.generate).pack(side='left',padx=6)
        ttk.Button(bar,text='Export CSV',command=self.export_samples).pack(side='left',padx=6)
        self.samples_label=ttk.Label(p,text='No samples generated');self.samples_label.pack(anchor='w',pady=12)
        tableframe=ttk.Frame(p);tableframe.pack(fill='both',expand=True)
        self.table=ttk.Treeview(tableframe,columns=('id','length','fragment'),show='headings')
        for k,label,width in [('id','Read #',80),('length','Length',90),('fragment','DNA fragment',750)]:
            self.table.heading(k,text=label);self.table.column(k,width=width,anchor='w')
        sc=ttk.Scrollbar(tableframe,orient='vertical',command=self.table.yview);self.table.configure(yscrollcommand=sc.set)
        self.table.pack(side='left',fill='both',expand=True);sc.pack(side='right',fill='y')

    def build_assembly(self):
        p=self.assembly_tab
        ttk.Label(p,text='Reconstruct the DNA',font=('Segoe UI',16,'bold')).pack(anchor='w')
        ttk.Label(p,text='Greedy exact suffix–prefix overlaps; no access to the original fragment positions.').pack(anchor='w',pady=(4,12))
        bar=ttk.Frame(p);bar.pack(fill='x')
        ttk.Label(bar,text='Minimum overlap (bases)').pack(side='left')
        ttk.Spinbox(bar,from_=5,to=110,textvariable=self.overlap_var,width=8).pack(side='left',padx=10)
        self.run_button=ttk.Button(bar,text='▶ Run reconstruction',command=self.run_assembly);self.run_button.pack(side='left',padx=8)
        ttk.Button(bar,text='Save report',command=self.save_report).pack(side='left',padx=8)
        self.results_label=ttk.Label(p,text='Generate samples, then run reconstruction.',font=('Segoe UI',11))
        self.results_label.pack(anchor='w',pady=15)
        self.results_text=tk.Text(p,height=15,wrap='word',font=('Consolas',11),bg='#111d2b',fg='#e1e8f1',padx=15,pady=15)
        self.results_text.pack(fill='both',expand=True)
        self.write_results('Assembly results will appear here.\n\nImportant limitations:\n• Repeated DNA can create misleading overlaps.\n• Random reads may leave gaps.\n• Greedy assembly may choose the wrong merge.\n• Comparing all pairs is computationally expensive.')

    def build_live(self):
        p=self.live_tab
        ttk.Label(p,text='Live reconstruction & replay',font=('Segoe UI',16,'bold')).pack(anchor='w')
        ttk.Label(p,text='Each frame shows an actual greedy merge. Recording runs in memory while assembly computes.').pack(anchor='w',pady=(4,12))
        controls=ttk.Frame(p);controls.pack(fill='x',pady=6)
        ttk.Button(controls,text='▶ Play',command=self.play_replay).pack(side='left',padx=3)
        ttk.Button(controls,text='⏸ Pause',command=self.pause_replay).pack(side='left',padx=3)
        ttk.Button(controls,text='◀ Previous',command=lambda:self.step_replay(-1)).pack(side='left',padx=3)
        ttk.Button(controls,text='Next ▶',command=lambda:self.step_replay(1)).pack(side='left',padx=3)
        ttk.Button(controls,text='Save recording (.txt)',command=self.export_recording).pack(side='left',padx=12)
        ttk.Label(controls,text='Delay (ms):').pack(side='left')
        self.delay=tk.IntVar(value=250)
        ttk.Combobox(controls,textvariable=self.delay,values=(50,100,250,500,1000),state='readonly',width=7).pack(side='left',padx=5)
        self.follow_live=tk.BooleanVar(value=True)
        ttk.Checkbutton(controls,text='Follow live merges',variable=self.follow_live).pack(side='left',padx=10)
        self.frame_label=ttk.Label(p,text='Run reconstruction to record merge events.')
        self.frame_label.pack(anchor='w',pady=8)
        self.timeline=ttk.Scale(p,from_=0,to=1,command=self.seek_replay)
        self.timeline.pack(fill='x',pady=(0,12))
        self.live_text=tk.Text(p,wrap='word',font=('Consolas',12),bg='#111d2b',fg='#e1e8f1',padx=16,pady=16)
        self.live_text.pack(fill='both',expand=True)
        self.live_text.tag_configure('overlap',foreground='#46d4a6')
        self.live_text.tag_configure('heading',foreground='#64a9ff')
        self.render_replay()

    def render_replay(self):
        self.live_text.configure(state='normal');self.live_text.delete('1.0','end')
        if not self.history:
            self.live_text.insert('end','No recorded merges yet.\n\nGenerate samples, then click Run reconstruction.\nThe recording will populate as merges are found.')
            self.frame_label.configure(text='0 merges recorded')
        else:
            i=min(self.play_index,len(self.history)-1)
            _,step,a,b,n,merged,remaining,elapsed=self.history[i]
            self.frame_label.configure(text=f'Merge {step:,} / {len(self.history):,} recorded  |  Overlap: {n} bases  |  Contigs left: {remaining:,}  |  {elapsed:.2f}s')
            def add(label,value):
                self.live_text.insert('end',label+'\n','heading')
                self.live_text.insert('end',value+'\n\n')
            self.live_text.insert('end','FRAGMENT A\n','heading')
            self.live_text.insert('end',a[:-n]);self.live_text.insert('end',a[-n:],'overlap');self.live_text.insert('end','\n\n')
            self.live_text.insert('end','FRAGMENT B\n','heading')
            self.live_text.insert('end',b[:n],'overlap');self.live_text.insert('end',b[n:]+'\n\n')
            add('MATCHING OVERLAP',a[-n:])
            add('MERGED CONTIG',merged)
            self.timeline.configure(to=max(1,len(self.history)-1))
            self.timeline.set(i)
        self.live_text.configure(state='disabled')

    def seek_replay(self,value):
        if not self.history:return
        i=max(0,min(int(float(value)),len(self.history)-1))
        if i!=self.play_index:
            self.follow_live.set(False)
            self.play_index=i;self.render_replay()

    def step_replay(self,delta):
        self.pause_replay()
        self.follow_live.set(False)
        if self.history:
            self.play_index=max(0,min(self.play_index+delta,len(self.history)-1))
            self.render_replay()

    def play_replay(self):
        if self.playing or not self.history:return
        if self.play_index>=len(self.history)-1:self.play_index=0
        self.playing=True;self.replay_tick()

    def replay_tick(self):
        if not self.playing:return
        self.render_replay()
        if self.play_index>=len(self.history)-1:
            self.playing=False;return
        self.play_index+=1
        self.play_after=self.root.after(max(20,self.delay.get()),self.replay_tick)

    def pause_replay(self):
        self.playing=False
        if self.play_after is not None:
            self.root.after_cancel(self.play_after);self.play_after=None

    def export_recording(self):
        if not self.history:
            messagebox.showinfo('No recording','Run reconstruction first.');return
        path=filedialog.asksaveasfilename(defaultextension='.txt',filetypes=[('Text recording','*.txt')])
        if not path:return
        with open(path,'w',encoding='utf-8') as f:
            for _,step,a,b,n,merged,remaining,elapsed in self.history:
                f.write(f'MERGE {step} | overlap={n} | remaining={remaining} | elapsed={elapsed:.2f}s\nA: {a}\nB: {b}\nOverlap: {a[-n:]}\nMerged: {merged}\n\n')
        messagebox.showinfo('Saved','Reconstruction recording saved.')

    def write_results(self,text):
        self.results_text.configure(state='normal');self.results_text.delete('1.0','end');self.results_text.insert('1.0',text);self.results_text.configure(state='disabled')
    def browse(self):
        path=filedialog.askopenfilename(filetypes=[('Genome files','*.gz *.fna *.fa *.fasta *.txt'),('All files','*.*')])
        if path:self.load(Path(path))
    def load(self,path):
        try:name,seq,n=read_fasta(path)
        except Exception as exc:messagebox.showerror('Cannot open genome',str(exc));return
        self.name=name;self.genome=seq;self.path=str(path)
        self.info.set(f'{Path(path).name}  |  {len(seq):,} bases  |  {n} FASTA record(s) (using first)')
        self.select_region(False)
    def select_region(self,whole):
        if not self.genome:return
        self.seq=self.genome if whole else self.genome[:2000]
        self.reads=[];self.contigs=[];self.history=[];self.play_index=0;self.pause_replay();self.render_replay()
        self.samples_label.set('No samples generated') if isinstance(self.samples_label,tk.StringVar) else self.samples_label.configure(text='No samples generated')
        for item in self.table.get_children():self.table.delete(item)
        self.region_label.configure(text=f'Analyzing {len(self.seq):,} bases  •  Window at position 1')
        self.start_var.set(1);self.slider.configure(to=max(1,len(self.seq)-20));self.draw_sequence()
        self.results_label.configure(text='Generate samples, then run reconstruction.')
    def draw_sequence(self):
        if not self.seq:return
        pos=max(0,min(int(self.start_var.get())-1,len(self.seq)-1));length=max(1,min(int(self.window_var.get()),400))
        fragment=self.seq[pos:pos+length]
        self.region_label.configure(text=f'Selected region: {len(self.seq):,} bases   |   Viewing bases {pos+1:,}–{pos+len(fragment):,}')
        self.dna_text.configure(state='normal');self.dna_text.delete('1.0','end')
        for i,ch in enumerate(fragment):
            self.dna_text.insert('end',ch,ch if ch in COLORS else 'N')
            if (i+1)%50==0:self.dna_text.insert('end','\n')
            elif (i+1)%10==0:self.dna_text.insert('end',' ')
        self.dna_text.configure(state='disabled')
    def generate(self):
        if not self.seq:messagebox.showwarning('Missing genome','Load a genome first');return
        try:count=self.count_var.get();low=self.low_var.get();high=self.high_var.get();seed=self.seed_var.get()
        except tk.TclError:messagebox.showerror('Invalid values','Use whole numbers');return
        if not(1<=count<=10000 and 1<=low<=high<=len(self.seq)):
            messagebox.showerror('Invalid settings','Require 1–10,000 reads and 1 ≤ minimum ≤ maximum ≤ selected genome length.');return
        self.reads=make_reads(self.seq,count,low,high,seed);self.contigs=[];self.history=[];self.play_index=0;self.pause_replay();self.render_replay()
        for item in self.table.get_children():self.table.delete(item)
        for r in self.reads:self.table.insert('','end',values=(r['id'],len(r['sequence']),r['sequence']))
        covered=[False]*len(self.seq)
        for r in self.reads:
            for i in range(r['start'],r['start']+len(r['sequence'])):covered[i]=True
        self.coverage=sum(covered)/len(covered)*100
        self.samples_label.configure(text=f'{len(self.reads):,} reads generated  |  {self.coverage:.2f}% of positions sampled  |  {len(set(r["sequence"] for r in self.reads)):,} unique fragments')
        self.status.set('Samples generated. Ready to reconstruct.');self.tabs.select(self.samples_tab)
    def run_assembly(self):
        if self.busy:return
        if not self.reads:messagebox.showwarning('No samples','Generate samples first');return
        minimum=self.overlap_var.get()
        if minimum<1:messagebox.showerror('Invalid overlap','Minimum overlap must be positive');return
        self.history=[];self.play_index=0;self.pause_replay();self.follow_live.set(True);self.render_replay()
        self.busy=True;self.run_button.configure(state='disabled');self.progress.start(10);self.status.set('Assembling reads…')
        reads=list(self.reads)
        def worker():
            try:self.messages.put(('done',assemble(reads,minimum,self.messages.put)))
            except Exception as exc:self.messages.put(('error',str(exc)))
        threading.Thread(target=worker,daemon=True).start()
    def poll(self):
        try:
            while True:
                msg=self.messages.get_nowait()
                if msg[0]=='merge':
                    self.history.append(msg)
                    if self.follow_live.get() and not self.playing:
                        self.play_index=len(self.history)-1
                        if len(self.history)==1 or len(self.history)%5==0:
                            self.render_replay()
                elif msg[0]=='progress':self.status.set(f'Assembling: {msg[1]} merges, {msg[2]} contigs remaining, {msg[3]:.1f}s')
                elif msg[0]=='done':self.finish(*msg[1])
                elif msg[0]=='error':
                    self.busy=False;self.progress.stop();self.run_button.configure(state='normal');messagebox.showerror('Assembly failed',msg[1])
        except queue.Empty:pass
        self.root.after(120,self.poll)
    def finish(self,contigs,merges,elapsed):
        self.busy=False;self.progress.stop();self.run_button.configure(state='normal')
        self.contigs=contigs
        if self.follow_live.get() and self.history:self.play_index=len(self.history)-1
        self.render_replay()
        longest=contigs[0] if contigs else ''
        exact=longest==self.seq
        # A circular genome may reconstruct with a different rotation; this check is linear and exact.
        contained=self.seq.find(longest)>=0
        self.report=(f'Genome: {self.name}\nSelected reference length: {len(self.seq):,} bases\n'
                     f'Input reads: {len(self.reads):,}\nMinimum overlap: {self.overlap_var.get()} bases\n'
                     f'Coverage of sampled reference: {self.coverage:.2f}%\n'
                     f'Greedy merges: {merges:,}\nRemaining contigs: {len(contigs):,}\n'
                     f'Longest contig: {len(longest):,} bases\nElapsed time: {elapsed:.2f} seconds\n'
                     f'Exact reconstruction: {"YES" if exact else "NO"}\n'
                     f'Longest contig occurs exactly in reference: {"YES" if contained else "NO"}\n\n'
                     'CONTIG LENGTHS\n'+', '.join(str(len(x)) for x in contigs[:50])+ '\n\n'
                     'WHY THIS CAN FAIL\n'
                     '1. Repeated regions create ambiguous overlaps and incorrect merges.\n'
                     '2. Uncovered reference positions cannot be reconstructed from reads.\n'
                     '3. Greedy decisions are local and may not yield the true genome.\n'
                     '4. Pairwise comparisons are slow as the number of reads grows.\n'
                     '5. This model uses exact overlaps, no sequencing errors or reverse complements.\n'
                     '6. phiX174 is circular; a linear comparison may miss a correct rotation.\n\n'
                     'FIRST CONTIG (up to 1000 bases)\n'+longest[:1000])
        self.write_results(self.report)
        self.results_label.configure(text=f'{len(contigs)} contig(s)  |  longest {len(longest):,} bases  |  exact match: {"YES" if exact else "NO"}')
        self.status.set(f'Finished in {elapsed:.2f}s');self.tabs.select(self.assembly_tab)
    def export_samples(self):
        if not self.reads:messagebox.showinfo('No samples','Generate samples first');return
        path=filedialog.asksaveasfilename(defaultextension='.csv',filetypes=[('CSV','*.csv')])
        if path:
            import csv
            with open(path,'w',newline='',encoding='utf-8') as f:
                writer=csv.writer(f);writer.writerow(['read_id','length','sequence'])
                for r in self.reads:writer.writerow([r['id'],len(r['sequence']),r['sequence']])
            messagebox.showinfo('Saved','Reads exported without their original positions.')
    def save_report(self):
        if not self.contigs:messagebox.showinfo('No results','Run reconstruction first');return
        path=filedialog.asksaveasfilename(defaultextension='.txt',filetypes=[('Text report','*.txt')])
        if path:Path(path).write_text(self.report+'\n\nALL CONTIGS\n'+'\n'.join(self.contigs),encoding='utf-8')

if __name__=='__main__':
    root=tk.Tk();App(root);root.mainloop()
