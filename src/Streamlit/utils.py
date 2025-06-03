from supabase import create_client, Client
import streamlit as st
import pandas as pd
SUPABASE_URL = "https://unwqzwvqjdcwmlwafxoo.supabase.co"
SUPABASE_KEY = "eyJhbGciOiJIUzI1NiIsInR5cCI6IkpXVCJ9.eyJpc3MiOiJzdXBhYmFzZSIsInJlZiI6InVud3F6d3ZxamRjd21sd2FmeG9vIiwicm9sZSI6ImFub24iLCJpYXQiOjE3NDY4NjU0NzUsImV4cCI6MjA2MjQ0MTQ3NX0.7P3a3cL3JPUidguddJo7HFgy25Ft60UKPpfo4aznqlo"

@st.cache_resource
def load_supabase():
    url = st.secrets["supabase"]["url"]
    key = st.secrets["supabase"]["key"]
    return create_client(url, key)

supabase: Client = load_supabase()

def show_loading_animation(message="Loading data"):
    st.markdown(f"""
    <div class="loading-container">
        <div class="loading-spinner"></div>
        <div class="loading-text">{message}<span class="loading-dots">...</span></div>
    </div>
    """, unsafe_allow_html=True)

def show_skeleton_loader():
    st.markdown("""
    <div class="dataframe-container">
        <div style="padding: 1rem;">
            <div class="skeleton-row"></div>
            <div class="skeleton-row"></div>
            <div class="skeleton-row"></div>
            <div class="skeleton-row"></div>
            <div class="skeleton-row"></div>
            <div class="skeleton-row"></div>
            <div class="skeleton-row"></div>
            <div class="skeleton-row"></div>
        </div>
    </div>
    """, unsafe_allow_html=True)

def fetch_page(offset: int, limit: int = 50):
    try:
        response = supabase.table("Nucleotide") \
            .select("*") \
            .range(offset, offset + limit - 1) \
            .execute()
        return pd.DataFrame(response.data)
    except Exception as e:
        st.error(f"Error fetching data: {str(e)}")
        return pd.DataFrame()

def get_total_count():
    """Get total count of records for pagination info - optimized to prevent timeout"""
    try:
        response = supabase.table("Nucleotide") \
            .select("id", count="exact") \
            .limit(1) \
            .execute()
        return response.count
    except Exception as e:
        st.warning("⚠️ Could not get exact count. Showing approximate pagination.")
        return None
