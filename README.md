# Resource Allocation for Text Semantic Communications

This repository contains MATLAB code for resource allocation optimization in text semantic communication systems using DeepSC (Deep Learning-based Semantic Communication).

## Overview

This project implements and simulates resource allocation strategies for semantic communication networks, comparing a proposed semantic-aware model with conventional communication systems. The implementation uses the Hungarian algorithm for optimal channel assignment and evaluates system performance under various network conditions.

## Features

- **Semantic Communication Simulation**: Implements DeepSC-based semantic communication with configurable parameters
- **Resource Allocation Optimization**: Uses the Hungarian algorithm for optimal channel-user assignment
- **Multi-Model Comparison**: Compares proposed semantic-aware model with conventional 4G/5G systems
- **Performance Visualization**: Generates figures for semantic similarity, spectral efficiency, and system performance
- **Channel Modeling**: Includes realistic path loss, shadow fading, and Rayleigh fading models

## Requirements

- MATLAB (tested on R2019a or later)
- No additional toolboxes required

## Files Description

### Main Simulation Scripts

- **`fig_3.m`**: Main simulation comparing proposed semantic-aware model with conventional models
  - Simulates semantic spectral efficiency (S-SE) vs. number of channels
  - Compares proposed model with conventional models using different semantic symbol numbers
  - Uses Monte Carlo simulation (1000 iterations) for statistical accuracy

- **`fig2.m`**: Visualizes the semantic similarity performance table for DeepSC
  - Creates 3D surface plot showing semantic similarity vs. SNR and semantic symbols
  - Useful for understanding DeepSC performance characteristics

- **`fig4a_channels.m`**: Compares different communication systems across varying channel numbers
  - Evaluates semantic communication, ideal, 4G, and 5G systems
  - Analyzes how spectral efficiency scales with available channels

- **`fig4b.m`**: Studies the impact of transmission power on system performance
  - Compares semantic-aware vs. conventional systems across power range (-40 to 30 dBm)
  - Demonstrates power efficiency of semantic communication

- **`fig4c_transform.m`**: Analyzes the effect of transformation factor on system performance
  - Studies how semantic-to-bit conversion factor impacts spectral efficiency

### Supporting Files

- **`Hungarian.m`**: Implementation of the Hungarian algorithm for optimal assignment
  - Finds minimum cost matching in bipartite graphs
  - Used for optimal channel-user allocation
  - Based on the classical Kuhn-Munkres algorithm

- **`sem_table.mat`**: Pre-computed DeepSC performance lookup table
  - Contains semantic similarity values for different SNR and symbol configurations
  - Rows: semantic symbol numbers (1-20)
  - Columns: SNR values (-10 to 20 dB)

- **`Resource Allocation for Text Semantic Communications(1).pdf`**: Research paper with theoretical background and detailed analysis

## System Parameters

### Common Parameters
- **Number of devices (users)**: 5
- **Cell radius**: 500 meters
- **Noise power**: -121 dBm (180000 × 10^-17.4 mW)
- **Shadow fading factor**: 6 dB
- **Transmission power**: 10 dBm (configurable in some scripts)
- **Monte Carlo iterations**: 1000
- **Semantic similarity threshold**: 0.9
- **Transformation factor**: 40

### Model-Specific Parameters
- **Proposed Model**: Variable semantic symbols (1-20)
- **Conventional Model**: Fixed semantic symbols [3, 5, 7, 9, 11]
- **Spectral efficiency threshold**: 1 bit/s/Hz (conventional systems)

## Usage

### Running the Main Comparison Simulation

```matlab
% Run the main simulation (Figure 3)
fig_3

% This will:
% - Simulate both proposed and conventional models
% - Generate S-SE vs. number of channels plot
% - Save results for analysis
```

### Visualizing DeepSC Performance

```matlab
% Visualize semantic similarity table (Figure 2)
fig2

% Creates a 3D surface plot showing:
% - X-axis: Semantic symbols
% - Y-axis: SNR (dB)
% - Z-axis: Semantic similarity
```

### Comparing Different Systems

```matlab
% Compare systems across different channel numbers (Figure 4a)
fig4a_channels

% Compare systems across different transmission powers (Figure 4b)
fig4b

% Analyze transformation factor impact (Figure 4c)
fig4c_transform
```

## How It Works

### 1. Channel Model
The system models wireless channels with three components:
- **Path Loss**: 128.1 + 37.6 × log10(d/1000), where d is distance in meters
- **Shadow Fading**: Log-normal distribution with 6 dB standard deviation
- **Small-Scale Fading**: Rayleigh fading (exponential distribution)

### 2. Semantic Communication
- Uses DeepSC lookup table (`sem_table.mat`) for semantic similarity vs. SNR
- Semantic Spectral Efficiency (S-SE) = Semantic Similarity / Number of Symbols
- Threshold-based quality assurance (similarity ≥ 0.9)

### 3. Resource Allocation
- **Proposed Model**: Optimizes semantic symbol count per user-channel pair
- **Optimization**: Hungarian algorithm finds optimal channel-user assignment
- **Objective**: Maximize total semantic spectral efficiency (S-SE)

### 4. Performance Metrics
- **S-SE (Semantic Spectral Efficiency)**: suts/s/Hz
- **Conventional SE**: bits/s/Hz (Shannon capacity)
- **Channel Assignment**: Binary allocation matrix (1 = assigned, 0 = not assigned)

## Results

The simulations generate several figures:

1. **Figure 2**: 3D visualization of DeepSC semantic similarity performance
2. **Figure 3**: S-SE comparison between proposed and conventional models vs. channel count
3. **Figure 4a**: Multi-system comparison (Semantic, Ideal, 4G, 5G) vs. channels
4. **Figure 4b**: System performance vs. transmission power
5. **Figure 4c**: Impact of transformation factor on performance

Key findings (typical results):
- Proposed semantic-aware model outperforms conventional fixed-symbol approaches
- Performance scales with number of available channels
- Semantic communication shows better efficiency at moderate SNR ranges
- Optimal resource allocation significantly improves overall system performance

## Algorithm Workflow

```
1. Initialize simulation parameters
2. For each Monte Carlo iteration:
   a. Generate random user positions in cell
   b. Calculate large-scale channel gains (path loss + shadow fading)
   c. Generate small-scale fading (Rayleigh)
   d. For each user-channel pair:
      - Calculate SNR
      - Look up semantic similarity from table
      - Compute S-SE for all possible symbol counts
      - Select optimal symbol count (proposed model)
   e. Run Hungarian algorithm for channel assignment
   f. Calculate total system S-SE
3. Average results over all iterations
4. Plot performance metrics
```

## Customization

### Modifying System Parameters

```matlab
% In any script, you can modify:
n_devices = 5;          % Change number of users
radius = 500;           % Change cell radius (meters)
p = 10;                 % Change transmission power (dBm)
mento = 1000;           % Change number of Monte Carlo iterations
f_th = 0.9;             % Change semantic similarity threshold
```

### Adding New Models

To add a new communication model:
1. Define the spectral efficiency calculation method
2. Implement channel assignment strategy
3. Add to comparison loop in main scripts
4. Update plotting section with new curve

## Performance Optimization

- **Reduce simulation time**: Decrease `mento` (Monte Carlo iterations)
- **Increase accuracy**: Increase `mento` (1000 is good balance)
- **Parallel processing**: MATLAB's Parallel Computing Toolbox can speed up simulations

## Citation

If you use this code in your research, please cite the accompanying paper:

```
Resource Allocation for Text Semantic Communications
[Add full citation details from the PDF]
```

## License

Please refer to the repository license or contact the authors for usage terms.

## Contributing

For questions, suggestions, or contributions, please open an issue or contact the repository maintainer.

## Acknowledgments

- DeepSC framework for semantic communication
- Hungarian algorithm implementation by Alex Melin (2006)

## References

1. DeepSC: Deep Learning-based Semantic Communication System
2. Hungarian Algorithm for Assignment Problems
3. Wireless Channel Models (Path Loss, Shadow Fading, Rayleigh Fading)

## Contact

For more information, please refer to the PDF document included in this repository or contact the authors.

---

**Note**: This code is for research and educational purposes. The simulation parameters are based on typical wireless communication scenarios and can be adjusted for specific use cases.
