import streamlit as st
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
def apply_custom_css():
    st.markdown(f"""
    <style>
    /* Import Google Fonts */
    @import url('https://fonts.googleapis.com/css2?family=Inter:wght@300;400;500;600;700&display=swap');
    
    /* Global font settings */
    * {{
        font-family: {THEME_CONFIG['font_family']} !important;
        color: {THEME_CONFIG['text_color']} !important;
    }}
    
    /* Streamlit specific overrides */
    .stApp {{
        background-color: {THEME_CONFIG['background_color']} !important;
        color: {THEME_CONFIG['text_color']} !important;
    }}
    
    /* Remove default Streamlit padding */
    .main {{
        padding: 0;
        background-color: {THEME_CONFIG['background_color']} !important;
    }}
    
    /* Full width body */
    .stApp {{
        margin: 0;
        padding: 0;
    }}
    
    /* Main container styling */
    .main .block-container {{
        width: {THEME_CONFIG['container_width']};
        max-width: {THEME_CONFIG['max_container_width']};
        margin: 0;
        padding: {THEME_CONFIG['padding']};
        background-color: {THEME_CONFIG['background_color']};
    }}
    
    /* Header styling */
    .main-header {{
        text-align: center;
        color: {THEME_CONFIG['text_color']} !important;
        font-family: {THEME_CONFIG['font_family']} !important;
        margin-bottom: 2rem;
        padding: 1rem;
        background-color: {THEME_CONFIG['secondary_background_color']};
        border-radius: {THEME_CONFIG['border_radius']};
    }}
    
    .page-info{{
        margin: 20px 20px;
    }}
    .main-header h1 {{
        color: {THEME_CONFIG['text_color']} !important;
        font-weight: 700 !important;
    }}
    
    .main-header p {{
        color: {THEME_CONFIG['text_color']} !important;
        opacity: 0.8;
    }}
    
    /* Dataframe container */
    .dataframe-container {{
        background-color: {THEME_CONFIG['secondary_background_color']};
        padding: 1rem;
        border-radius: {THEME_CONFIG['border_radius']};
        margin: 1rem 0;
        box-shadow: {THEME_CONFIG['box_shadow']};
    }}
    
    /* Dataframe styling */
    .stDataFrame {{
        background-color: {THEME_CONFIG['secondary_background_color']} !important;
    }}
    
    .stDataFrame table {{
        background-color: {THEME_CONFIG['secondary_background_color']} !important;
        color: {THEME_CONFIG['text_color']} !important;
    }}
    
    .stDataFrame th {{
        background-color: {THEME_CONFIG['primary_color']} !important;
        color: white !important;
        font-weight: 600 !important;
    }}
    
    .stDataFrame td {{
        background-color: {THEME_CONFIG['secondary_background_color']} !important;
        color: {THEME_CONFIG['text_color']} !important;
    }}
    
    /* Info box styling */
    .stAlert {{
        background-color: {THEME_CONFIG['secondary_background_color']} !important;
        color: {THEME_CONFIG['text_color']} !important;
        border: none !important;
    }}
    
    .stAlert > div {{
        color: {THEME_CONFIG['text_color']} !important;
    }}
    
    /* Fix info background color */
    div[data-testid="stAlert"] {{
        background-color: {THEME_CONFIG['secondary_background_color']} !important;
        border: none !important;
    }}
    
    div[data-testid="stAlert"] div {{
        background-color: {THEME_CONFIG['secondary_background_color']} !important;
        color: {THEME_CONFIG['text_color']} !important;
        border: none !important;
    }}
    
    .nav-container {{
        display: flex;
        flex-wrap: wrap;
        justify-content: space-between;
        align-items: center;
        margin-top: 2rem;
        padding: 1rem;
        gap: 1rem;
        background-color: {THEME_CONFIG['secondary_background_color']};
        border-radius: {THEME_CONFIG['border_radius']};
    }}

    /* Responsive layout for small screens */
    @media (max-width: 768px) {{
        .nav-container {{
            flex-direction: column;
            gap: 1rem;
        }}
        
        /* Add spacing between button containers on mobile */
        .stButton {{
            margin: 0.5rem 0 !important;
            width: 100% !important;
        }}
        
        /* Ensure page info has proper spacing */
        .page-info {{
            margin: 1rem 0 !important;
            padding: 1rem !important;
        }}
    }}

    .stButton > button {{
        background-color: {THEME_CONFIG['primary_color']} !important;
        color: white !important;
        border: none !important;
        border-radius: {THEME_CONFIG['border_radius']} !important;
        padding: 0.5rem 1rem !important;
        font-family: {THEME_CONFIG['font_family']} !important;
        font-weight: 500 !important;
        transition: all 0.3s ease !important;
        width: 100% !important;
        min-height: 48px !important; /* Better touch target */
    }}

    .stButton > button:hover {{
        background-color: {THEME_CONFIG['primary_color']}dd !important;
        transform: translateY(-2px) !important;
        box-shadow: 0 4px 8px rgba(0, 0, 0, 0.2) !important;
    }}

    /* Mobile-specific button styling */
    @media (max-width: 768px) {{
        .stButton > button {{
            margin-bottom: 0.75rem !important;
            min-height: 52px !important; /* Larger touch target on mobile */
            font-size: 1rem !important;
        }}
    }}
    
    /* Page info styling */
    .page-info {{
        text-align: center;
        color: {THEME_CONFIG['text_color']} !important;
        font-family: {THEME_CONFIG['font_family']} !important;
        font-weight: 500;
        background-color: {THEME_CONFIG['secondary_background_color']};
        padding: 0.5rem 1rem;
        border-radius: {THEME_CONFIG['border_radius']};
        margin: 0 1rem;
    }}
    
    /* Loading Animation Styles */
    .loading-container {{
        display: flex;
        flex-direction: column;
        align-items: center;
        justify-content: center;
        padding: 3rem;
        background-color: {THEME_CONFIG['secondary_background_color']};
        border-radius: {THEME_CONFIG['border_radius']};
        margin: 1rem 0;
    }}
    
    .loading-spinner {{
        width: 50px;
        height: 50px;
        border: 4px solid {THEME_CONFIG['secondary_background_color']};
        border-top: 4px solid {THEME_CONFIG['primary_color']};
        border-radius: 50%;
        animation: spin 1s linear infinite;
        margin-bottom: 1rem;
    }}
    
    @keyframes spin {{
        0% {{ transform: rotate(0deg); }}
        100% {{ transform: rotate(360deg); }}
    }}
    
    .loading-text {{
        color: {THEME_CONFIG['text_color']} !important;
        font-family: {THEME_CONFIG['font_family']} !important;
        font-size: 1.1rem;
        font-weight: 500;
        opacity: 0.8;
    }}
    
    .loading-dots {{
        display: inline-block;
        animation: dots 1.5s infinite;
    }}
    
    @keyframes dots {{
        0%, 20% {{
            color: rgba(255, 255, 255, 0);
            text-shadow:
                .25em 0 0 rgba(255, 255, 255, 0),
                .5em 0 0 rgba(255, 255, 255, 0);
        }}
        40% {{
            color: {THEME_CONFIG['text_color']};
            text-shadow:
                .25em 0 0 rgba(255, 255, 255, 0),
                .5em 0 0 rgba(255, 255, 255, 0);
        }}
        60% {{
            text-shadow:
                .25em 0 0 {THEME_CONFIG['text_color']},
                .5em 0 0 rgba(255, 255, 255, 0);
        }}
        80%, 100% {{
            text-shadow:
                .25em 0 0 {THEME_CONFIG['text_color']},
                .5em 0 0 {THEME_CONFIG['text_color']};
        }}
    }}
    
    /* Pulse animation for skeleton loading */
    .skeleton-row {{
        background: linear-gradient(90deg, {THEME_CONFIG['secondary_background_color']} 25%, #3A3A3A 50%, {THEME_CONFIG['secondary_background_color']} 75%);
        background-size: 200% 100%;
        animation: pulse 1.5s infinite;
        height: 40px;
        margin: 8px 0;
        border-radius: 4px;
    }}
    
    @keyframes pulse {{
        0% {{
            background-position: 200% 0%;
        }}
        100% {{
            background-position: -200% 0%;
        }}
    }}
    .stButton[data-testid="baseButton-secondary"] > button {{
        background: linear-gradient(45deg, #57268C, #7B4FBF) !important;
        color: white !important;
        border: 2px solid #57268C !important;
        border-radius: 0.5rem !important;
        padding: 0.7rem 1.5rem !important;
        font-weight: 600 !important;
        transition: all 0.3s ease !important;
        box-shadow: 0 2px 4px rgba(87, 38, 140, 0.3) !important;
    }}

    .stButton[data-testid="baseButton-secondary"] > button:hover {{
        background: linear-gradient(45deg, #7B4FBF, #57268C) !important;
        transform: translateY(-2px) !important;
        box-shadow: 0 4px 8px rgba(87, 38, 140, 0.4) !important;
        border-color: #7B4FBF !important;
    }}

    /* Dialog styling */
    div[data-testid="stModal"] {{
        background-color: rgba(26, 26, 26, 0.95) !important;
    }}

    div[data-testid="stModal"] > div {{
        background-color: #2A2A2A !important;
        border: 1px solid #57268C !important;
        border-radius: 1rem !important;
        box-shadow: 0 10px 25px rgba(0, 0, 0, 0.5) !important;
        color: #FFFFFF !important;
    }}

    /* File uploader styling */
    .stFileUploader > div {{
        background-color: #2A2A2A !important;
        border: 2px dashed #57268C !important;
        border-radius: 0.5rem !important;
        color: #FFFFFF !important;
    }}

    .stFileUploader label {{
        color: #FFFFFF !important;
    }}

    /* Expander styling */
    .streamlit-expanderHeader {{
        background-color: #2A2A2A !important;
        color: #FFFFFF !important;
    }}

    .streamlit-expanderContent {{
        background-color: #1A1A1A !important;
        color: #FFFFFF !important;
    }}
    /* Hide Streamlit menu and footer */
    #MainMenu {{visibility: hidden;}}
    footer {{visibility: hidden;}}
    header {{visibility: hidden;}}
    
    /* Custom scrollbar */
    ::-webkit-scrollbar {{
        width: 8px;
    }}
    ::-webkit-scrollbar-track {{
        background: {THEME_CONFIG['secondary_background_color']};
    }}
    ::-webkit-scrollbar-thumb {{
        background: {THEME_CONFIG['primary_color']};
        border-radius: 4px;
    }}
    </style>
    """, unsafe_allow_html=True)