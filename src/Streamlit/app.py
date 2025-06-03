import streamlit as st
import pandas as pd
from css import apply_custom_css
from  utils import fetch_page, get_total_count, show_loading_animation, show_skeleton_loader
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
            nav_col1, nav_col2, nav_col3 = st.columns([2, 1, 2])
            
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
            
            st.markdown('</div>', unsafe_allow_html=True)
        else:
            if st.session_state.page > 0:
                st.button("Go Back", key="back_btn", on_click=lambda: handle_page_change(-1))

def handle_page_change(direction):
    """Handle page navigation with loading state"""
    st.session_state.page += direction
    st.session_state.loading = True

if __name__ == "__main__":
    main()