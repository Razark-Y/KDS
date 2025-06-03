# KDS: Trait-Based SARS-CoV-2 Lineage Classification

This repository contains the source code and analysis pipeline for classifying SARS-CoV-2 lineages (Alpha, Delta, Omicron) using raw FASTA sequences—without relying on curated metadata or alignment-based mutation detection.

---

## 🔗 Live Demo

Access the Streamlit web interface here:  
👉 [https://nucleotide-prediction.streamlit.app/](https://nucleotide-prediction.streamlit.app/)

You can input your FASTA sequence and receive a predicted lineage in real time.

---

## 📁 Repository Structure

- `data/` — Contains the curated dataset used for training and validation.  
- `Prediction/` — Contains the inference logic and trained model pipeline.  
- `Analysis/` — Includes evaluation scripts, visualizations, and feature analysis.  
- `Format/` — Responsible for sequence scraping and formatting scripts.  
- `Preprocessing/` — Contains aligned FASTA files and related preprocessing tools using the Wuhan reference genome.
- `Streamlit/` — Contains interface for inference.

---

## 🔐 Supabase Access

This project connects to a Supabase database using a **static URL and API key**.  
⚠️ **Note:** The key may change if the current instance exceeds its egress limit.  
If connection issues occur, please update the key or check the latest configuration.

---

## 📄 Citation

If you use this project or its ideas, please cite:

**"Trait-Based Differentiation of SARS-CoV-2 Variants: A PANGOLIN-Inspired Approach to Analyzing Alpha, Delta, and Omicron"**

---

Maintained by  
👨‍💻 [Wilson Yusda](mailto:13522019@std.stei.itb.ac.id) | Bandung Institute of Technology
👨‍💻 [Enrique Yanuar](mailto:13522077@std.stei.itb.ac.id) | Bandung Institute of Technology
👨‍💻 [Mesach Harmasendro](mailto:135220117@std.stei.itb.ac.id) | Bandung Institute of Technology

