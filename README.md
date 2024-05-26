## Study, Implementation, and Simulation of OFDM and FBMC Systems in MATLAB
# Overview
This project explores and implements two key multicarrier modulation schemes: Orthogonal Frequency Division Multiplexing (OFDM) and Filter Bank Multicarrier (FBMC) systems, specifically using Offset Quadrature Amplitude Modulation (OQAM). Both systems are simulated and analyzed using MATLAB to evaluate their performance, particularly in the context of 5G communication systems.

Table of Contents
Introduction to Modulation Schemes for 5G
1.1 Basics of OFDM
1.2 Further Multicarrier Modulation Schemes
1.2.1 Filter Bank Multi-Carrier (FBMC)
1.2.2 Generalized Frequency Division Multiplexing (GFDM)
1.3 Objectives
Detailed Analysis of OFDM
2.1 Advantages
2.2 Disadvantages
2.3 Implementation in MATLAB
2.4 Simulation Results
Detailed Analysis of FBMC
3.1 FBMC Physical Layer
3.2 Implementation using PHYDYAS Filter in MATLAB
3.2.1 Prototype Filter Design
3.2.2 OQAM
3.2.3 Frequency Spreading
3.2.4 Polyphase Network Implementation
3.3 Simulation Results
3.4 Comparison between OFDM and OQAM-FBMC
Conclusions
Introduction to Modulation Schemes for 5G
Multicarrier modulation is crucial for 5G communications. OFDM has been the dominant technology for broadband multicarrier systems. However, for applications like cognitive radios and uplink of multiuser systems, OFDM might not be ideal. FBMC is considered a promising alternative due to its efficient use of spectrum and reduced intercarrier interference (ICI).

Basics of OFDM
OFDM splits a high-rate data stream into several lower rate streams transmitted over multiple subcarriers, reducing delay spread and eliminating intersymbol interference (ISI) through guard intervals. Key parameters include the number of subcarriers, guard time, symbol duration, and subcarrier spacing.

Further Multicarrier Modulation Schemes
Filter Bank Multi-Carrier (FBMC)
FBMC uses subcarrier filtering with OQAM, providing better spectral efficiency and reduced out-of-band emissions compared to OFDM. It does not require a cyclic prefix, saving bandwidth and power.

Generalized Frequency Division Multiplexing (GFDM)
GFDM is a flexible multicarrier method without strict orthogonality, suitable for cognitive radio applications. It manages out-of-band emissions well but has a more complex receiver design.

Objectives
The project aims to:

Understand OFDM and FBMC systems.
Implement and simulate these systems in MATLAB.
Compare their performance, particularly in terms of Bit Error Rate (BER) and Power Spectral Density (PSD).
Detailed Analysis of OFDM
Advantages
High spectral efficiency.
Robustness against frequency-selective fading.
Simplified channel equalization.
Disadvantages
High Peak-to-Average Power Ratio (PAPR).
Sensitivity to frequency offset and phase noise.
Implementation in MATLAB
The OFDM system is implemented using MATLAB, focusing on the IFFT and FFT blocks for modulation and demodulation, respectively. Guard intervals are added to mitigate ISI.

Simulation Results
Simulations involve evaluating BER and PSD across different SNR levels. The results show expected performance aligning with theoretical predictions.

Detailed Analysis of FBMC
FBMC Physical Layer
FBMC uses subcarrier filtering and OQAM to improve spectral efficiency and reduce out-of-band emissions.

Implementation using PHYDYAS Filter in MATLAB
Prototype Filter Design
Designing the prototype filter involves setting the filter length and analyzing its frequency response.

OQAM
OQAM splits complex symbols into real and imaginary parts, transmitted with a time offset to maintain orthogonality.

Frequency Spreading
Frequency spreading improves spectral containment and reduces interference.

Polyphase Network Implementation
The polyphase network efficiently implements the filter bank, reducing computational complexity compared to the IFFT approach in OFDM.

Simulation Results
Simulations show FBMC's performance in terms of BER and PSD, highlighting its advantages in spectral efficiency and reduced interference.

Comparison between OFDM and OQAM-FBMC
FBMC demonstrates lower sidelobe levels and better spectral efficiency without the need for a cyclic prefix, making it more suitable for certain 5G applications.

Conclusions
The study concludes that while OFDM is well-suited for many broadband applications, FBMC offers significant advantages in spectral efficiency and interference management, particularly for uplink multiuser and cognitive radio scenarios.

This README provides a high-level summary of the detailed study, implementation, and simulation of OFDM and FBMC systems in MATLAB. For in-depth analysis and MATLAB code, please refer to the full document.
