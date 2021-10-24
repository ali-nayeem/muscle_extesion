# Multiobjective Application-aware MUSCLE
We develop Multiobjective Application-aware MUSCLE by embedding the following four application-aware objective functions, 
identified based on their better correlation to the tree accuracy, within the iterative phase of MUSCLE:
1. Maximize similarity for columns containing gaps (SIMG)
2. Maximize similarity for columns containing no gaps (SIMNG)
3. Maximize sum-of-pairs (SP)
4. Minimize number of gaps (GAP)

## Multiobjective principles
The goal of an multiobjective (MO) algorithm is to generate a set of solutions, popularly known as the Pareto-optimal solutions in the MO literature, 
which represent the best compromise among the (conflicting) objectives. Among the several classes of MO algorithms (e.g., pareto-based, decomposition-based, 
indicator-based, etc.), decomposition-based strategies are found effective to face the difficulties in handling 'many' (i.e., more than three) objectives. 
These algorithms decompose the task of generating several alternative solutions into many single-objective problems with the help of a set of well-distributed 
weight vectors, popularly known as reference directions. Each weight vector aggregates the different objective scores into a single value that eventually 
leads to one member of the final solution set.

## Simplified workflow
We drive the iterative search process of MUSCLE with a total of four objectives directed by a 4D weight vector. Figure 1 depicts a high-level workflow for one weight vector, where the steps (3.4 to 3.6) inspired by the MO approach are marked as red. This workflow is executed for all weight vectors to get alternative solutions and can be performed independently in parallel.
Also note that, unlike original MUSCLE, we update the guide tree (step 3.6) each time a better MSA is obtained to intensify the effect of MO principles.

<figure>
  <img src="https://github.com/ali-nayeem/mammle/blob/macos/diagram/ma-muscle.png" alt="Trulli" style="width:80%">
  <figcaption>Fig.1: High-level workflow of multiobjective application-aware MUSCLE for a single weight vector. Steps (3.4 to 3.6) added/modified on the original MUSCLE are marked with red color. This figure is modification of the original image taken from https://doi.org/10.1093/nar/gkh340.</figcaption>
</figure>

## Installation 
The current version has been developed and tested on Linux and MAC. 

You need to have:

- [GNU C++ Compiler](https://gcc.gnu.org/) 

Open a terminal and clone our [github repository](https://github.com/ali-nayeem/muscle_extesion). For example you can use: 

```bash
git clone https://github.com/ali-nayeem/muscle-extension.git
cd muscle_extesion
make
```   
## Execution

The weight values for the four objective functions can be passed via the following options:

```bash
./muscle -simg <w1> -simng <w2> -osp <w3> -gap <w4> -objscore sp -in input_seq.fasta -out out.fasta
```
