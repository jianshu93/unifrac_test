
C.F lucblassel/phylotree#11 and lucblassel/phylotree#8

# UniFrac implememtation in rust

This is an example repo to show how to compute the [UniFrac]() distance between a pair of samples containing taxa. 
It uses the [`phylotree`]() crate to parse the tree file and compute the unifrac distance.

Given a tree with 4 tips *(in file `tree_rot.nwk`)*, this program computes the unifrac distance between samples: 
 - $S_A = \{T_1,T_2,T_3,T_4\}$
 - $S_B = \{T_2,T_3\}$

Program output: 
```
 ❯ cargo run --release -- test_rot.nwk 
   Compiling unifrac v0.1.0 (/Users/lucblassel/Development/rust/unifrac)
    Finished `release` profile [optimized] target(s) in 0.71s
     Running `target/release/unifrac test_rot.nwk`
Tree: ((T1:0.2,(T2:0.1,T3:0.4):0.3):0.5,T4:0.6);
Leaf order = ["T1", "T2", "T3", "T4"]
B =
[[1, 1, 1, 1],
 [1, 1, 1, 0],
 [1, 0, 0, 0],
 [0, 1, 1, 0],
 [0, 1, 0, 0],
 [0, 0, 1, 0],
 [0, 0, 0, 1]]
b = [0, 0.5, 0.2, 0.3, 0.1, 0.4, 0.6]

P_A:[1, 1, 1, 1, 1, 1, 1]
P_B:[1, 1, 0, 1, 1, 1, 0]

Computing unifrac:
P_A.T @ Br @ P_B   = 1.3        (49ns)
Sum(P_A * P_B * b) = 1.3        (32ns)
D_unifrac(S_A,S_B) = 0.381
```
