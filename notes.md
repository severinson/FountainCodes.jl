
The goal is to run decoding experiments for Raptor and LDPC codes
However, I know that the current implementation is faulty, since it sometimes will report a decoding failure for full rank matrices
Further, I don't like the current API
I'd prefer to have it be a solver in the style of A \ b for finite field matrices
First, I need to make sure that my decoder is correct
Let's change the API a bit
The user has a SparseMatrixCSC, where each column corresponds to a constraint
From that matrix, the user creates a Decoder object, which stores the state of the decoder
Next, the decode operation passes around both the sparse matrix and the decoder object

There are some questions regarding decoding, since Raptor codes are over finite fields, whereas we want to use floats
It may be the case that using inactivation decoding is misguided, since there are efficient iterative least squares solvers,
e.g., Kazmarz (?)
But inactivation decoding is the route that they've chosen

Erasure code implementations and the decoder are too tangled up
The erasure code should just create generator (or parity check) matrices
The decoder should just solve linear systems
I can fix that if I start by re-writing the decoder to just solve linear systems


TODO
- Remove the phase and status from the decoder object
- Rename the metric names
- Check store_metrics everywhere needed
- Don't need to store num_symbols separately, since it's encoded by the matrix size
- Add Ti as a type parameter, and use that type consistently