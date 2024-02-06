CC = mpicc -Wall
CFLAGS = -O 
CLIBS = -lm

INCS = ac2dr_aux.h ac2dr_config.h
PROG = ac2dr

SRC = ac2dr_config.c \
		ac2dr_main.c \
	  ac2dr_aux.c \
	  ac2dr_input.c \
	  ac2dr_RK4.c \
	  ac2dr_RK4_damping.c \
	  ac2dr_BC.c \
	  ac2dr_getSBP6.c \
	  ac2dr_output.c \
	  ac2dr_src.c \
	  ac2dr_init.c \
	  ac2dr_dqdt2_sym_sg_ghost.c \
	  ac2dr_dqdt4_sym_sg_ghost.c \
	  ac2dr_sg.c \
	  ac2dr_material.c \
	  ac2dr_decompose_dom.c \
	  ac2dr_mpi_aux.c \
	  ac2dr_dqdt12.c \
	  ac2dr_dqdt6.c \
	  ac2dr_dqdt6_damping.c \
	  ac2dr_dqdt4.c 

OBJ = ${SRC:.c=.o}

.c.o:
	${CC} ${CFLAGS} -c $<

all: ${PROG}

${PROG}: ${OBJ}
	${CC} ${CFLAGS} -o ${PROG} ${OBJ} ${CLIBS}

clean:
	rm *.o ${PROG}
