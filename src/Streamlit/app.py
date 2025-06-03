import streamlit as st
import pandas as pd
from css import apply_custom_css
from utils import fetch_page, get_total_count, show_loading_animation, show_skeleton_loader
from preprocessor import *
import sys
import types

def setup_custom_classes():
    if '__main__' not in sys.modules:
        sys.modules['__main__'] = types.ModuleType('__main__')
    from custom_transformers import NRemover, TwoMerPCA, ThreeMerPCA, FeatureCreator, FeatureDropper
    print("Setting up custom classes in __main__")
    sys.modules['__main__'].NRemover = NRemover
    sys.modules['__main__'].TwoMerPCA = TwoMerPCA
    sys.modules['__main__'].ThreeMerPCA = ThreeMerPCA
    sys.modules['__main__'].FeatureCreator = FeatureCreator
    sys.modules['__main__'].FeatureDropper = FeatureDropper

THEME_CONFIG = {
    "primary_color": "#57268C",
    "background_color": "#1A1A1A",
    "secondary_background_color": "#2A2A2A",
    "text_color": "#FFFFFF",
    "font_family": "-apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif",
    "border_radius": "0.5rem",
    "container_width": "100vw",
    "max_container_width": "100%",
    "padding": "2rem",
    "box_shadow": "0 4px 6px -1px rgba(0, 0, 0, 0.1)"
}

@st.dialog("🔬 Perform Prediction")
def prediction_dialog():
    st.markdown("""
    <div style="text-align: center; margin-bottom: 1.5rem;">
        <h3 style="color: #57268C; margin-bottom: 0.5rem;">AI Prediction Analysis</h3>
        <p style="opacity: 0.8;">Enter the details below to perform nucleotide sequence prediction</p>
    </div>
    """, unsafe_allow_html=True)
    
    col1, col2 = st.columns(2)
    
    with col1:
        accession_id = st.text_input(
            "🧬 Accession ID",
            placeholder="Enter accession ID (e.g., NC_000001)",
            help="Unique identifier for the nucleotide sequence"
        )
    
    with col2:
        location = st.text_input(
            "📍 Location", 
            placeholder="Enter location (e.g., USA, China)",
            help="Geographic location where the sample was collected"
        )
    
    st.markdown("---")
    
    st.markdown("**📄 FASTA File Upload**")
    uploaded_file = st.file_uploader(
        "Choose a FASTA file",
        type=['txt', 'fasta', 'fa', 'fas'],
        help="Upload a text file containing FASTA sequence data",
        label_visibility="collapsed"
    )
    
    if uploaded_file is not None:
        file_details = {
            "Filename": uploaded_file.name,
            "File size": f"{uploaded_file.size} bytes",
            "File type": uploaded_file.type if uploaded_file.type else "text/plain"
        }
        st.success("✅ File uploaded successfully!")
        with st.expander("📋 File Details"):
            for key, value in file_details.items():
                st.write(f"**{key}:** {value}")
        with st.expander("👀 Preview File Content"):
            try:
                file_content = uploaded_file.read().decode('utf-8')
                preview_content = file_content[:500]
                if len(file_content) > 500:
                    preview_content += "..."
                st.code(preview_content, language="text")
                st.caption(f"Showing first 500 characters of {len(file_content)} total characters")
                uploaded_file.seek(0)
            except Exception as e:
                st.error(f"Error reading file: {str(e)}")
    st.markdown("---")
    col1, col2, col3 = st.columns([1, 2, 1])
    with col2:
        if st.button("🚀 Run Prediction", type="primary", use_container_width=True):
            if not accession_id.strip():
                st.error("❌ Please enter an Accession ID")
                return
            
            if not location.strip():
                st.error("❌ Please enter a Location")
                return
                
            if uploaded_file is None:
                st.error("❌ Please upload a FASTA file")
                return

            with st.spinner("🔄 Running AI prediction..."):
                import tempfile
                with tempfile.NamedTemporaryFile(delete=False, suffix=".fasta") as tmp:
                    tmp.write(uploaded_file.read())
                    temp_path = tmp.name
                setup_custom_classes()
                prediction_class = formatting_and_prediction(
                    file_path=temp_path,
                    location=location,
                    accession_id=accession_id,
                    collection_date="2020-01-01" 
                )

                label_map = {0: "Alpha", 1: "Delta", 2: "Omicron"}
                predicted_label = label_map.get(prediction_class, "Unknown")

                st.session_state.prediction_result = {
                    "accession_id": accession_id,
                    "location": location,
                    "filename": uploaded_file.name,
                    "file_size": uploaded_file.size,
                    "status": "completed",
                    "prediction": predicted_label
                }

            st.success("✅ Prediction completed successfully!")
            st.balloons()
            st.markdown("### 📊 Prediction Results")
            result_data = {
                "Accession ID": accession_id,
                "Location": location,
                "File": uploaded_file.name,
                "Prediction": predicted_label,
                "Status": "✅ Completed"
            }

            for key, value in result_data.items():
                st.write(f"**{key}:** {value}")

def main():
    apply_custom_css()
    
    if "page" not in st.session_state:
        st.session_state.page = 0
    if "loading" not in st.session_state:
        st.session_state.loading = False
    
    st.markdown("""
    <div class="main-header">
        <h1>🧬 Nucleotide Records Viewer</h1>
        <p>Browse nucleotide sequence data with advanced pagination</p>
    </div>
    """, unsafe_allow_html=True)
    
    data_placeholder = st.empty()
    nav_placeholder = st.empty()
    
    if "current_data" not in st.session_state or st.session_state.loading:
        with data_placeholder.container():
            show_loading_animation("Fetching nucleotide records")
            show_skeleton_loader()
        
        offset = st.session_state.page * 50
        df = fetch_page(offset)
        total_count = get_total_count()
        
        st.session_state.current_data = df
        st.session_state.total_count = total_count
        st.session_state.loading = False
        
        data_placeholder.empty()
    
    df = st.session_state.current_data
    total_count = st.session_state.get('total_count', None)
    offset = st.session_state.page * 50
    
    selected_columns = ["accession_id", "collection_date", "location", "status"]
    
    with data_placeholder.container():
        if not df.empty:
            start_record = offset + 1
            end_record = min(offset + len(df), offset + len(df))
            if total_count is not None:
                st.info(f"📊 Showing records {start_record:,}-{end_record:,} of {total_count:,}")
            else:
                st.info(f"📊 Showing records {start_record:,}-{end_record:,}")
            
            st.dataframe(
                df[selected_columns],
                use_container_width=True,
                hide_index=True
            )
        else:
            st.markdown("""
            <div class="dataframe-container">
                <div style="text-align: center; padding: 2rem;">
                    <h3>📭 No Data Available</h3>
                    <p>No records found for the current page.</p>
                </div>
            </div>
            """, unsafe_allow_html=True)
    
    with nav_placeholder.container():
        if not df.empty:
            st.markdown('<div style="margin-top: 1rem;">', unsafe_allow_html=True)
            
            nav_col1, nav_col2, nav_col3 = st.columns([1, 3, 1])
            
            with nav_col1:
                if st.session_state.page > 0:
                    st.button("Previous", key="prev_btn", on_click=lambda: handle_page_change(-1))
                else:
                    st.markdown("""
                    <div style="
                        background-color: rgba(87, 38, 140, 0.3);
                        color: white;
                        text-align: center;
                        padding: 0.7rem 0;
                        border-radius: 0.5rem;
                        opacity: 0.4;
                        cursor: not-allowed;
                    ">
                        Previous
                    </div>
                    """, unsafe_allow_html=True)
            
            with nav_col2:
                current_page = st.session_state.page + 1
                if total_count is not None:
                    total_pages = (total_count + 49) // 50 
                    st.markdown(f"""
                    <div class="page-info">
                        Page {current_page} of {total_pages} • {total_count:,} total records
                    </div>
                    """, unsafe_allow_html=True)
                else:
                    st.markdown(f"""
                    <div class="page-info">
                        Page {current_page} • Showing {len(df)} records
                    </div>
                    """, unsafe_allow_html=True)
            
            with nav_col3:
                if len(df) == 50:  
                    st.button("Next", key="next_btn", on_click=lambda: handle_page_change(1))
            
            st.markdown('<div style="margin-top: 1rem;">', unsafe_allow_html=True)
            
            pred_col1, pred_col2, pred_col3 = st.columns([1, 1, 1])
            with pred_col2:
                if st.button("🔬 Perform Prediction", key="pred_btn", type="secondary", use_container_width=True):
                    prediction_dialog()
            
            st.markdown('</div>', unsafe_allow_html=True)
            st.markdown('</div>', unsafe_allow_html=True)
        else:
            if st.session_state.page > 0:
                st.button("Go Back", key="back_btn", on_click=lambda: handle_page_change(-1))

def handle_page_change(direction):
    st.session_state.page += direction
    st.session_state.loading = True

if __name__ == "__main__":
    main()
