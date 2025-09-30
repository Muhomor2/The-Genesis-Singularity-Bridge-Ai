# Ω³ Protocol – Invisible Watermarking for AI Models

The **Ω³ Protocol** is a mathematically rigorous and empirically validated method for embedding **invisible watermarks** into the floating-point weights of neural networks.

Unlike metadata-based approaches that can be easily removed, Ω³ modifies the mantissa bits of selected parameters. This creates a **distributed, stealth signature** that remains hidden within the model’s numerical representation.

---

## 🔑 Key Features

* **Accuracy Preservation**
  Embedding does not affect model accuracy (ΔTop-1 ≈ –0.01% on ResNet-18).

* **Mathematical Rigor**
  Perturbation is strictly bounded:
  [
  \Delta x \leq 2^{-4} \cdot ULP(x),
  ]
  which equals ≈10⁻⁷ for FP32 and ≈10⁻¹⁵ for FP64.

* **High Stealth**
  Detection probability:
  [
  \Psi_{stealth} \approx 4 \times 10^{-15}.
  ]

* **Robustness**
  Survives local noise injection and pruning.
  Vulnerabilities remain under quantization (FP32 → INT8) and format conversion (e.g., JSON export).

---

## 📖 Documentation

* **`AI_MODEL_WATERMARK_OMEGA3.md`** – Full academic report:

  * Floating-point steganography background
  * Embedding operator and mathematical bounds
  * Empirical validation on ResNet-18 (ImageNet-1k)
  * Attack and robustness analysis
  * Comparison with baseline (metadata hash)
  * Future research directions

---

## 🚀 Applications

* **Intellectual property protection** for deep learning models
* **Provenance verification** of distributed AI artifacts
* **Federated learning watermarking**
* **Potential extension** to quantum amplitude encoding

---

## 📌 Citation

If you reference this work, please cite:

**Chechelnitsky, I. (2025). "The Ω³ Protocol: Mathematical Foundation, Empirical Validation, and Robustness Analysis."**

---
