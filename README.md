# LDPC Decoder in AWGN Channel Simulation
>This project provides a comprehensive simulation environment for (1023, 781) LDPC (Low-Density Parity-Check) decoding in the presence of an Additive White Gaussian Noise (AWGN) channel using the Sum-Product Algorithm (SPA).
>
## Simulation Flow
>**1.Encoding Phase:** The original message undergoes LDPC encoding.
>
>**2.Transmission Phase:** The encoded sequence gets transmitted via an AWGN channel, where it might be subject to noise and potential data corruption.
>
>**3.Decoding Phase:** The decoder aims to recover the original message from the received sequence, proceeding through several stages:
>
>>* Bottom-up phase: The check nodes process the incoming messages.
>>
>>* Top-down phase: The variable nodes handle the incoming messages.
>>
>>* Termination phase: This phase concludes with the final decisions being made for the estimated codeword.
>>
## Optimization in Decoding
>**Tree Reduction Algorithm:** This technique is applied in the bottom-up phase to combine messages efficiently.
>
>**Minimizing Redundant Computations:** During the decoding process, a significant amount of redundant computations can arise when updating different nodes.<br>
To address this, the decoder first calculates the overall result without excluding the node that needs to be updated. Once this result is obtained, the impact of the node to be updated is removed from the total, and then the cleaned-up value is finally updated.<br>
This strategy not only eliminates redundant computations but also speeds up the decoding process.
>
## Files in this Repository:
>**LDPC.c:** The main source code of the simulation.
>
>**ldpc_G_1023.txt:** Represents the generator matrix essential for the (1023, 781) LDPC code.
<div align="center">
<img src=https://github.com/wzaq4350/LDPC-code/blob/main/ldpc_G_1023.png width=60% />
</div>

>**ldpc_H_1023.txt:** This is a 2046x32 matrix as shown below. Rows 1-1023 record the connections of check nodes 1-1023 to their respective 32 variable nodes. Rows 1024-2046 document the links of variable nodes 1-1023 to their corresponding 32 check nodes.
>
<div align="center">
<img src=https://github.com/wzaq4350/LDPC-code/blob/main/ldpc_H_1023.jpg width=60% />
</div>
