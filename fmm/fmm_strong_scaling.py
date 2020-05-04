from __future__ import print_function, division

import json
import ctypes

import getopt
import sys
import math
import ppmd as md
import ppmd.cuda as mdc
import time
from ppmd.coulomb.fmm import *

from scipy import constants

kB = constants.Boltzmann



READ = md.access.READ
WRITE = md.access.WRITE
INC_ZERO = md.access.INC_ZERO
INC = md.access.INC

rank =  md.mpi.MPI.COMM_WORLD.Get_rank()
nproc = md.mpi.MPI.COMM_WORLD.Get_size()



shell_steps = 8
shell_width = 0.25
ewald_cutoff = 30
atom_types = ('Na', 'Cl')

DISABLE_PRINTING = False
EWALD = False
FMM = True
L = 10

scale_factor = 0.00010364


PositionDat = md.data.PositionDat
ParticleDat = md.data.ParticleDat
ScalarArray = md.data.ScalarArray
State = md.state.State
ParticleLoop = md.loop.ParticleLoop
Pairloop = md.pairloop.PairLoopNeighbourListNSOMP
#Pairloop = md.pairloop.CellByCellOMP
PBC = md.domain.BoundaryTypePeriodic()
ScalarArray = md.data.ScalarArray
GlobalArray = md.data.GlobalArray

prefix = os.path.abspath(os.path.expanduser(sys.argv[1]))
parameters = json.loads(open(sys.argv[2]).read())

CONFIG  = os.path.join(prefix, 'CONFIG')
CONTROL = os.path.join(prefix, 'CONTROL')
FIELD   = os.path.join(prefix, 'FIELD')

extent = md.utility.dl_poly.read_domain_extent(CONFIG)

if rank == 0:
    positions = np.load(os.path.join(prefix, 'positions.npy'))
    velocities = np.load(os.path.join(prefix, 'velocities.npy'))
    charges = np.load(os.path.join(prefix, 'charges.npy'))



get_con = md.utility.dl_poly.get_control_value
get_fie = md.utility.dl_poly.get_field_value

dt = float(get_con(CONTROL, 'TIMESTEP')[0][0])
dt /= math.sqrt(scale_factor)

# steps = int(get_con(CONTROL, 'STEPS')[0][0])
steps = parameters['steps']

cutoff = float(get_con(CONTROL, 'RVDW')[0][0])
ewald_precision = float(get_con(CONTROL, 'EWALD PRECISION')[0][0])

N = sum([int(nx[0]) for nx in get_fie(FIELD, 'NUMMOLS')])

type0 = get_fie(FIELD, atom_types[0])
mass0 = float(type0[0][0])
charge0 = float(type0[0][1])
type1 = get_fie(FIELD, atom_types[1])
mass1 = float(type1[0][0])
charge1 = float(type1[0][1])


atoms_types_map ={
    atom_types[0]: 0,
    atom_types[1]: 1
}

mass_map = {
    atom_types[0]: mass0,
    atom_types[1]: mass1
}

charge_map = {
    atom_types[0]: charge0,
    atom_types[1]: charge1
}


lj00 = get_fie(FIELD, '{} {} LJ'.format(atom_types[0], atom_types[0]))
lj00_eps = float(lj00[0][0])
lj00_sigma = float(lj00[0][1])
lj11 = get_fie(FIELD, '{} {} LJ'.format(atom_types[1], atom_types[1]))
lj11_eps = float(lj11[0][0])
lj11_sigma = float(lj11[0][1])
lj01 = get_fie(FIELD, '{} {} LJ'.format(atom_types[0], atom_types[1]))
lj01_eps = float(lj01[0][0])
lj01_sigma = float(lj01[0][1])

# ----------------------------------------------------------------------
if rank == 0:
    print(60*'-')
    print('extent', extent)
    print('N', N)
    print('dt', dt)
    print('steps', steps)
    print('shell_steps', shell_steps)
    print('cutoff', cutoff)
    print('ewald_precision', ewald_precision)
    print('t0', atom_types[0])
    print('m0', mass0)
    print('q0', charge0)
    print('t1', atom_types[1])
    print('m1', mass1)
    print('q1', charge1)
    print('0-0', lj00_eps, lj00_sigma)
    print('1-1', lj11_eps, lj11_sigma)
    print('0-1', lj01_eps, lj01_sigma)
    print(60*'-')
# ----------------------------------------------------------------------


# create a state called A
A = State()
A.npart = N
rn = cutoff + shell_width

# Init the domain in A
A.domain = md.domain.BaseDomainHalo(extent=extent)
A.domain.boundary_condition = PBC


# init particle dats for md part
A.p = PositionDat(ncomp=3)
A.p0 = ParticleDat(ncomp=3)

A.v = ParticleDat(ncomp=3)
A.f = ParticleDat(ncomp=3)
A.mass = ParticleDat(ncomp=1)
A.q = ParticleDat(ncomp=1)
A.type = ParticleDat(ncomp=1, dtype=ctypes.c_int)

A.u   = GlobalArray(ncomp=1)
A.ke  = GlobalArray(ncomp=1)

if EWALD:
    A.crr = ScalarArray(ncomp=1)
    A.cri = ScalarArray(ncomp=1)
    A.crs = ScalarArray(ncomp=1)


if rank == 0:
    A.p[:] = positions[:]
    velocities *= math.sqrt(scale_factor)
    


    for px in range(N):
        charge = charges[px]
        if charge < 0:
            sym = 'Cl'
        else:
            sym = 'Na'

        A.type[px, 0] = atoms_types_map[sym]
        A.mass[px, 0] = mass_map[sym]
        A.q[px, 0] = charge

    A.v[:] = velocities[:]


A.f[:] = 0.0
A.u.set(0.0)
A.ke.set(0.0)
if EWALD:
    A.crr[0] = 0.0
    A.cri[0] = 0.0


A.scatter_data_from(0)

# ---------------------------------------
# create a pairwise kernel and pairloop

# 2 species kernel that should vectorise assuming types in {0, 1}
lj_kernel_code = '''
const double R0 = P.j[0] - P.i[0];
const double R1 = P.j[1] - P.i[1];
const double R2 = P.j[2] - P.i[2];
const double r2 = R0*R0 + R1*R1 + R2*R2;

const double dt0 = (double) T.i[0];
const double dt1 = (double) T.j[0];
const double dt0dt1 = dt0*dt1;

// epsilon
double e = e00;
e += dt0 * e01me00;
e += dt1 * e01me00;
e += dt0dt1 * e_coeff;

// sigma^2
double s = s00;
s += dt0 * s01ms00;
s += dt1 * s01ms00;
s += dt0dt1 * s_coeff;

// potential shift
double ljs = ljs00;
ljs += dt0 * ljs01mljs00;
ljs += dt1 * ljs01mljs00;
ljs += dt0dt1 * ljs_coeff;

// avoid the divide
double eos = eos00;
eos += dt0 * eos01meos00;
eos += dt1 * eos01meos00;
eos += dt0dt1 * eos_coeff;

// compute the interaction
const double r_m2 = s/r2;
const double r_m4 = r_m2*r_m2;
const double r_m6 = r_m4*r_m2;

u[0] += (r2 < rc2) ? 2. * e * ((r_m6-1.0)*r_m6 + ljs) : 0.0;

const double r_m8 = r_m4*r_m4;
const double f_tmp = (-48.0 * eos) * (r_m6 - 0.5) * r_m8;

F.i[0]+= (r2 < rc2) ? f_tmp*R0 : 0.0;
F.i[1]+= (r2 < rc2) ? f_tmp*R1 : 0.0;
F.i[2]+= (r2 < rc2) ? f_tmp*R2 : 0.0;
'''

const = md.kernel.Constant

shift00 = (lj00_sigma/cutoff)**6. - (lj00_sigma/cutoff)**12.
shift01 = (lj01_sigma/cutoff)**6. - (lj01_sigma/cutoff)**12.
shift11 = (lj11_sigma/cutoff)**6. - (lj11_sigma/cutoff)**12.

# shift00 = 0.
# shift01 = 0.
# shift11 = 0.


eoscoeff = -2.0 * lj01_eps/(lj01_sigma**2.0) + lj00_eps/(lj00_sigma**2.0) + \
    lj11_eps/(lj11_sigma**2.0)

constants = (
    const('e00', lj00_eps),
    const('e01me00', lj01_eps - lj00_eps),
    const('e_coeff', -2.*(lj01_eps - lj00_eps) + lj11_eps - lj00_eps),

    const('s00', lj00_sigma**2.0),
    const('s01ms00', lj01_sigma**2.0 - lj00_sigma**2.0),
    const('s_coeff', -2.*(lj01_sigma**2.0 - lj00_sigma**2.0) + \
          lj11_sigma**2.0 - lj00_sigma**2.0),

    const('eos00', lj00_eps/(lj00_sigma**2.0)),
    const('eos01meos00', lj01_eps/(lj01_sigma**2.0) - lj00_eps/(lj00_sigma**2.0)),
    const('eos_coeff', eoscoeff),

    const('ljs00', shift00),
    const('ljs01mljs00', shift01 - shift00),
    const('ljs_coeff', -2.*(shift01 - shift00) + shift11 - shift00),

    const('rc2', cutoff**2.)
)

ljkernel = md.kernel.Kernel('two_species_lj', lj_kernel_code, constants)

# pairloop to update forces and potential energy
vdw_force = Pairloop(
    kernel=ljkernel,
    dat_dict={
        'P':A.p(md.access.READ),
        'T':A.type(md.access.READ),
        'F':A.f(md.access.INC_ZERO),
        'u':A.u(md.access.INC_ZERO)
    },
    shell_cutoff=cutoff+shell_width
)

# ---------------------------------------
# create a velocity verlet particle loop
# kernel and create a particle loop object

vv_kernel1_code = '''
const double M_tmp = dht/M.i[0];
V.i[0] += F.i[0]*M_tmp;
V.i[1] += F.i[1]*M_tmp;
V.i[2] += F.i[2]*M_tmp;
P.i[0] += dt*V.i[0];
P.i[1] += dt*V.i[1];
P.i[2] += dt*V.i[2];
'''

vv_kernel2_code = '''
//printf("%f %f %f\\n", V.i[0], V.i[1], V.i[2]);
const double M_tmp = dht/M.i[0];
V.i[0] += F.i[0]*M_tmp;
V.i[1] += F.i[1]*M_tmp;
V.i[2] += F.i[2]*M_tmp;
k[0] += (V.i[0]*V.i[0] + V.i[1]*V.i[1] + V.i[2]*V.i[2])*0.5*M.i[0];
'''


constants = [
    const('dt', dt),
    const('dht', 0.5*dt),
    const('N', N)
]

vv_kernel1 = md.kernel.Kernel('vv1', vv_kernel1_code, constants)
vv_p1 = ParticleLoop(
    kernel=vv_kernel1,
    dat_dict={'P': A.p(md.access.W),
              'V': A.v(md.access.W),
              'F': A.f(md.access.R),
              'M': A.mass(md.access.R)}
)

vv_kernel2 = md.kernel.Kernel('vv2', vv_kernel2_code, constants)
vv_p2 = ParticleLoop(
    kernel=vv_kernel2,
    dat_dict={'V': A.v(md.access.W),
              'F': A.f(md.access.R),
              'M': A.mass(md.access.R),
              'k': A.ke(md.access.INC0)}
)




# ---------------------------------------
# ewald class instance

if EWALD:

    c = md.coulomb.ewald.EwaldOrthoganal(
        domain=A.domain,
        eps=ewald_precision,
        real_cutoff=ewald_cutoff*0.5,
        shared_memory='omp',
        force_unit=internal_to_ev(),
        energy_unit=internal_to_ev()
    )

    c.evaluate_contributions(positions=A.p, charges=A.q)
    A.cri[0] = 0.0
    c.extract_forces_energy_reciprocal(A.p, A.q, A.f, A.cri)
    A.crr[0] = 0.0
    c.extract_forces_energy_real(A.p, A.q, A.f, A.crr)
    A.crs[0] = 0.0
    c.evaluate_self_interactions(A.q, A.crs)




if FMM:
    #R = int(max(log(0.8*N,10), 3))
    R = 6
    fmm = PyFMM(domain=A.domain, r=R, l=L, eps=None, free_space=False,
            force_unit=internal_to_ev(), energy_unit=internal_to_ev())

    if rank == 0:
        print("internal to ev", internal_to_ev())
        print("R =", R)
# ---------------------------------------
# use the integrator range class to implement
# a timestepping loop


t0 = time.time()
ke_list = []
u_list = []
cou_list = []


# convert to eV
qfactor = 1.0
if EWALD:
    qfactor = c.internal_to_ev()
dltfactor = 0.001957528385665669

tfactor = (2.0 * 1.602176565e-19) /(3. * N * kB)

if rank == 0:
    print("kelvin factor", tfactor)


def print_header():
    if rank == 0:
        print("{:^8s} {:^16s} {:^16s} {:^16s} {:^16s} {:^16s} {:^10s}".format(
            "it "," KE (K)"," U (eV)"," Cou (eV)", "cfg (eV)"," total (eV)", "time (s)"))


ewald_pot = 0.0
fmm_pot = 0.0

dump_index = 0


print_every = max(steps - 2, 1)

for it in md.method.IntegratorRange(steps, dt, A.v, shell_steps, shell_width,
                                    verbose=True):

    # velocity verlet 1
    vv_p1.execute(A.npart_local)
    
    A.f[:N:, 0] = 0
    vdw_force.execute()

    # ewald
    if EWALD:
        c.evaluate_contributions(positions=A.p, charges=A.q)
        A.cri[0] = 0.0
        c.extract_forces_energy_reciprocal(A.p, A.q, A.f, A.cri)
        A.crr[0] = 0.0
        c.extract_forces_energy_real(A.p, A.q, A.f, A.crr)
        A.crs[0] = 0.0
        c.evaluate_self_interactions(A.q, A.crs)
        ewald_pot = A.crr[0] + A.cri[0] + A.crs[0]
    
    if FMM:
        fmm_pot = fmm(positions=A.p, charges=A.q, forces=A.f)
    
    #A.f[:N:, 0] = 0

    # velocity verlet 2
    vv_p2.execute(A.npart_local)


    if (it % print_every == 0) and not DISABLE_PRINTING:
        if FMM:
            cou_list.append(fmm_pot)
        elif EWALD:
            cou_list.append(ewald_pot)
        else:
            cou_list.append(0)

        ke_list.append(A.ke[0,])
        u_list.append(A.u[0])
        
        toteng = ke_list[-1] + u_list[-1]
        if FMM or EWALD:
            toteng += cou_list[-1]

        if rank == 0:
            if it % (print_every*60) == 0:
                print_header()
            print("{:8d} {:16.8e} {:16.8e} {:16.8e} {:16.8e} {:16.8e} {: 10.4f}".format(
                it,
                float(ke_list[-1]*tfactor),
                float(u_list[-1]),
                float(cou_list[-1]),
                float(cou_list[-1] + u_list[-1]),
                float( toteng ),
                time.time() - t0
            ))

total_time = time.time() - t0

if rank == 0:
    d = {
        'num_mpi_ranks': nproc,
        'num_omp_threads': runtime.NUM_THREADS,
        'time_taken': total_time,
        'time_taken_per_step': total_time / steps,
    }
    with open('last_time.json', 'w') as fh:
        fh.write(json.dumps(d))







