### Lab 2: Oliker--Prussner Method

#### Ting Lin, Peking University

```
- do.m					demo file
- loadfunction.m		load the function
- getequation.m			get the discrete setting
- meshinit.m			initial the mesh with Delaunay triangulation
- OPinit.m				Lift the boundary
- OPsolve.m				Perron's iteration
- quadpts.m				get the numerical integration (from iFEM)
- W2perror.m			discrete W2p error
- report/				report dir
```



The demo file:

```matlab
[u,f] = loadfunction(1);
h = 1/4;
[U,F, G, X,Y,bd, elem, elemind, adj, count, bdadj] = getequation(u,f,h);
[U, X,Y, bd, elem,elemind, adj, count] = OPinit(U,F, G, X,Y,bd, elem, elemind, adj, count,1000,bdadj, 1e-6);
[U, X,Y, bd, elem,elemind, adj, count]  = OPsolve(U,F, G,X,Y,bd,elem, elemind, adj, count, 10000,1e-6);
```

