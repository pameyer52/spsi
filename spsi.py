#!/usr/bin/env python

''' 
simplest possible structure interpolator 
'''

#see if we're in pymol or not
try:
    from pymol import cmd
    import os
    #we are, so setup a wrapper
    def morph(start_obj, end_obj, morph_obj, nsteps = 10 ):
        ''' pymol wrapper for spsi
        start_obj = object name for morph start
        end_obj = object name for morph end
        morph_obj = output object name
        nsteps = number of interpolation steps
        '''
        #setup tmp input files
        pid = os.getpid() #process id for unique temporary files names
        tmpfile_st = 'tmp1-%d.pdb' % pid
        tmpfile_en = 'tmp2-%d.pdb' % pid
        tmpfile_op = 'tmp3-%d.pdb' % pid
        cmd.do('save %s, %s' %(tmpfile_st, start_obj, ) )
        cmd.do('save %s, %s' %(tmpfile_en, end_obj, ) )

        #run
        spsi(tmpfile_st, tmpfile_en, tmpfile_op, nsteps)
        
        #load results
        cmd.do('load %s, %s' %(tmpfile_op, morph_obj, ) )

        #cleanup 
        os.remove(tmpfile_st)
        os.remove(tmpfile_en)
        os.remove(tmpfile_op)

except ImportError:
    #not in pymol, proceed as usual
    pass


import sqlite3

# support functions 
def parse_pdb_atom_line(ln):
    ''' parse a PDB formatted ATOM or HETATM line , return as dictionary '''
    if not ( ln.startswith('ATOM') or ln.startswith('HETATM') ):
        raise Exception('Not valid ATOM line')
    d = {}
    d['serial'] = int( ln[6:11] ) #atom serial number
    d['name'] = ln[12:16] #atom name
    d['resName'] = ln[17:20] #residue name
    d['chainID'] = ln[21] 
    d['resSeq'] = int( ln[22:26] ) #residue number
    d['x'] = float( ln[30:38] )
    d['y'] = float( ln[38:46] )
    d['z'] = float( ln[46:54] )
    d['q'] = float( ln[54:60] ) #occupancy
    d['b'] = float( ln[60:66] ) #b-factor
    return d

def generate_pdb_atom_line(serial, name, resName, chainID, resSeq, x, y, z, q = 1.0 , b = 23.42):
    ''' generate PDB formatted ATOM line '''
    fmt = 'ATOM %6d %s %s %s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n' 
    return fmt % (serial, name, resName, chainID, resSeq, x, y, z, q, b, )

def is_atom_line(ln):
    ''' true if line starts with ATOM or HETATM
    '''
    ia = ln.startswith('ATOM')
    ih = ln.startswith('HETATM')
    if ia or ih:
        return True
    return False

def load_pdb_atoms(pdbfile):
    ''' read ATOM or HETATM lines from a pdb file, return list of
    atom dictionaries
    '''
    inp = open(pdbfile,'r')
    atoms = [parse_pdb_atom_line(ln) for ln in inp.readlines() if \
    is_atom_line(ln) ]
    inp.close()
    return atoms

def store_atoms(db_cur, atom_lst, istep):
    ''' store atoms in db '''
    for atom in atom_lst:
        db_cur.execute('INSERT INTO atoms(atom_name, residue_number, chain_id, x, y, z, i_step, residue_name) values(?, ?, ?, ?,?,?, ?, ?)', ( atom['name'], atom['resSeq'], atom['chainID'], atom['x'], atom['y'], atom['z'], istep, atom['resName'],) )

def linear_interpolate(st, en, p):
    ''' simple linear interpolation between two coordinates by percentage '''
    d = en - st
    return st + (d * p )

#main routine 
def spsi(start_file, end_file, out_file, nsteps = 10):
    '''
    start_file = starting structure
    end_file = ending structure
    out_file = output (multi-model) pdb file
    nsteps = number of intermediate steps 
    '''
    print('loading input pdb files')
    #load both input files, store in db
    conn = sqlite3.connect(':memory:') #in memory DB, don't need persistence
    cur = conn.cursor()
    def sqlexec(c, sql):
        ''' wrapper for executing raw sql statement'''
        c.execute(sql)
    sqlexec(cur, 'CREATE TABLE atoms(id INTEGER PRIMARY KEY AUTOINCREMENT, atom_name TEXT,residue_number INTEGER,residue_name TEXT,chain_id TEXT,x REAL,y REAL,z REAL,i_step INTEGER,unpaired_flag INTEGER DEFAULT 1)')
    sqlexec(cur,'CREATE INDEX idx_atom_chid ON atoms(chain_id)')
    sqlexec(cur,'CREATE INDEX idx_atom_resn ON atoms(residue_number)')
    sqlexec(cur,'CREATE INDEX idx_atom_atnam ON atoms(atom_name)')
    sqlexec(cur,'CREATE TABLE pairs(id_start INTEGER,id_stop INTEGER)')
    start_model = load_pdb_atoms(start_file)
    end_model = load_pdb_atoms(end_file)
    store_atoms(cur, start_model, 0)
    store_atoms(cur, end_model, nsteps-1) 

    print('pairing atoms')
    #pair atoms by atom name, residue number and chain id.

    #views look like the simplest way of dealing with correct ordering of
    # start/end pairs.
    cur.execute('CREATE VIEW a AS SELECT id,residue_number, atom_name, chain_id FROM atoms WHERE i_step = 0')
    vstmt = 'CREATE VIEW b AS SELECT id,residue_number, atom_name, chain_id FROM atoms WHERE i_step = %d' % (nsteps-1, ) #if nsteps wasn't constrained to be integer, this would be a sql injection vunerability
    cur.execute(vstmt)
    cur.execute('SELECT a.id, b.id FROM a INNER JOIN b WHERE a.residue_number=b.residue_number AND a.atom_name=b.atom_name AND a.chain_id = b.chain_id and a.id <> b.id')
    r = cur.fetchall()
    for row in r:
        id_st = row[0]
        id_en = row[1]
        cur.execute('INSERT INTO pairs(id_start, id_stop) VALUES(?,?)',(id_st,id_en,) )
        cur.execute('UPDATE atoms SET unpaired_flag=0 WHERE id=?',(id_st,))
        cur.execute('UPDATE atoms SET unpaired_flag=0 WHERE id=?',(id_en,))

    #deal with unpaired atoms 
    cur.execute('DELETE FROM atoms WHERE unpaired_flag=1')

    print('interpolating coordinates')
    #create intermediate models by linear interpolation of coordinates
    cur.execute('SELECT id_start, id_stop FROM pairs')
    r = cur.fetchall()
    for row in r:
        id_st = row[0]
        id_en = row[1]
        cur.execute('SELECT atom_name, residue_number, chain_id, x, y, z,residue_name FROM atoms WHERE id=?',(id_st,))
        ar = cur.fetchall()[0] #should be only one result
        name = ar[0]
        resn = ar[1]
        chid = ar[2]
        x_st = ar[3]
        y_st = ar[4]
        z_st = ar[5]
        resname = ar[6]
        cur.execute('SELECT x,y,z FROM atoms WHERE id=?',(id_en,))
        ar2 = cur.fetchall()[0]
        x_en = ar2[0]
        y_en = ar2[1]
        z_en = ar2[2]
        #loop over intermediate structures
        for istep in range(1, nsteps-1):
            pct = float(istep) / nsteps
            #linear interpolation of x,y,z 
            x_intp = linear_interpolate( x_st, x_en, pct )
            y_intp = linear_interpolate( y_st, y_en, pct )
            z_intp = linear_interpolate( z_st, z_en, pct )
            cur.execute('INSERT INTO atoms(atom_name, residue_number, chain_id, x,y,z, i_step, residue_name) values(?,?,?, ?,?,?, ?, ?)',(name, resn, chid, x_intp,y_intp,z_intp, istep, resname, ) )

    print('writing output')
    #output multi-model pdb file
    opf = open(out_file, 'w')
    iatom = 0
    for imodel in range(0, nsteps):
        opf.write('MODEL  %7d\n'%(imodel,))
        cur.execute('SELECT atom_name,residue_name,chain_id,residue_number,x,y,z FROM atoms WHERE i_step=? ORDER BY chain_id, residue_number',(imodel,))
        for row in cur.fetchall():
            s = generate_pdb_atom_line(iatom, row[0], row[1], row[2], row[3], row[4], row[5], row[6], )
            opf.write(s)
            iatom += 1
        opf.write('ENDMDL\n')
    opf.close()
    cur.close()
    conn.commit()
    conn.close()
    print('done')


if __name__ == '__main__':
    import sys
    try:
        startfile = sys.argv[1]
        stopfile = sys.argv[2]
        morphfile = sys.argv[3]
        ns = int(sys.argv[4])
    except IndexError:
        print('format is: spsi.py [start pdb file] [end pdb file] [output (morph) pdb file] [number of interpolation steps]')
        sys.exit(1)
    spsi(startfile, stopfile, morphfile, ns)
