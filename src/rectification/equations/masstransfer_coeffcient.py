from rectification.utils import unitcheck
from scipy.constants import g
import numpy as np



@unitcheck(sigma_lc_top="N/m", sigma_hc_top="N/m", x_aver_top_mass="kmol/kmol", res_unit="N/m")
def sigma_top(sigma_lc_top, sigma_hc_top, x_aver_top_mass):
    """
    Calculates the surface tension at the top of column.
    Parameters
    ----------
    sigma_lc_top : float
        The surface tension of low-boilling component at the top of column, [N / m]
    sigma_hc_top : float
        The surface tension of high-boilling component at the top of column, [N / m]
    x_aver_top_mass : float
        The average mass concentration at top of column, [kg/kg]
    Returns
    -------
    sigma_top : float
        The surface tension at the top of column, [N / m]
    References
    ----------
    &&&&&
    """       
    return (sigma_lc_top * x_aver_top_mass  + (1 - x_aver_top_mass) * sigma_hc_top)


@unitcheck(sigma_lc_bot="N/m", sigma_hc_bot="N/m", x_aver_bot_mass="kmol/kmol", res_unit="N/m")
def sigma_bot(sigma_lc_bot, sigma_hc_bot, x_aver_bot_mass):
    """
    Calculates the surface tension at the bottom of column.
    Parameters
    ----------
    sigma_lc_bot : float
        The surface tension of low-boilling component at the bottom of column, [N / m]
    sigma_hc_bot : float
        The surface tension of high-boilling component at the bottom of column, [N / m]
    x_aver_bot_mass : float
        The average mass concentration at bot of column, [kg/kg]
    Returns
    -------
    sigma_bot : float
        The surface tension at the  bottom of column, [N / m]
    References
    ----------
    &&&&&
    """       
    return (sigma_lc_bot * x_aver_bot_mass  + (1 - x_aver_bot_mass) * sigma_hc_bot)


#region Calculating of bubble layer
@unitcheck(q_liq_top="m**2/s", h_septum="m", w_oper="m/s", mu_liq_top="Pa/s", sigma_mix="N/m", sigma_water="N/m", res_unit="m")
def heigth_layer_top(q_liq_top, h_septum, w_oper, m_coef, mu_liq_top, sigma_top, sigma_water_top):
    """
    Calculates the heigth ligth layer of  the liquid.
    Parameters
    ----------
    q_liq_top : float
        The specific flow rate of the liquid for 1 m of drain septum at the top of column, [m**2/s]
    h_septum : float
        The heigth of drain septum, [m]
    w_oper : float
        The operating speed in the column, [m/s]
    mu_liq_top : float
        The viscocity of mix liquid at the top of column, [Pa/s]
    m_coef : float
        The specific coefficient for this equation, [dimensionless]
    sigma_top : float
        The surface tension at the top of column, [N / m]
    sigma_water_top : float
        The surface tension of water [N/m]
    Returns
    -------
    heigth_layer_top : float
        The heigth ligth layer of  the liquid at the top of coulmn, [m]
    References
    ----------
    Дытнерский, страница 239, формула 6.39
    """
    return 0.787 * q_liq_top**0.2 * h_septum**0.56 * w_oper**(m_coef) * (1 - 0.31 * np.exp(-0.11 * mu_liq_top)) * (sigma_top/sigma_water_top)**0.09


@unitcheck(q_liq_bot="m**2/s", h_septum="m", w_oper="m/s", mu_liq_bot="Pa/s", sigma_mix="N/m", sigma_water="N/m", res_unit="m")
def heigth_layer_bot(q_liq_bot, h_septum, w_oper, m_coef, mu_liq_bot, sigma_bot, sigma_water_bot):
    """
    Calculates the heigth ligth layer of  the liquid at the bottom of column.
    Parameters
    ----------
    q_liq_bot : float
        The specific flow rate of the liquid for 1 m of drain septum at the bottoms of column, [m**2/s]
    h_septum : float
        The heigth of drain septum, [m]
    w_oper : float
        The operating speed in the column at the bottoms of column, [m/s]
    mu_liq_bot : float
        The viscocity of mix liquid at the bottoms of column[Pa/s]
    m_coef : float
        The specific coefficient for this equation, [dimensionless]
    sigma_bot : float
        The surface tension of mix at the bottoms of column, [N/m]
    sigma_water_bot : float
        The surface tension of water at the bottoms of column, [N/m]
    Returns
    -------
    heigth_layer : float
        The heigth ligth layer of  the liquid at the bottom of column, [m]
    References
    ----------
    Дытнерский, страница 239, формула 6.39
    """
    return 0.787 * q_liq_bot**0.2 * h_septum**0.56 * w_oper**m_coef * (1 - 0.31 * np.exp(-0.11 * mu_liq_bot)) * (sigma_bot/sigma_water_bot)**0.09


def m_coef_top(h_septum):
    """
    Calculates the specific coefficient for calculation the heigth ligth layer of liquid equation.
    Parameters
    ----------
    h_septum : float
        The heigth of drain septum, [m]
    Returns
    -------
    m_coef : float
        The specific coefficient for this equation [dimensionless]
    References
    ----------
    Дытнерский, страница 239, формула 6.39
    """
    return 0.05 - h_septum*4.6


@unitcheck(rho_top_liq="kg/m**3", L_top="kg/s", L_septum="m", res_unit="m**2/s")
def q_liq_top(rho_top_liq, L_septum, L_top):
    """
    Calculates the specific flow rate of the liquid for 1 m of drain septum
    Parameters
    ----------
    rho_top_liq : float
        The destiny of liquid at top of column, [kg/m**3]
    L_top : float
        The flow rate of liquid at top of column, [kg/s]
    L_septum : float
        The length of drain septum, [m]
    Returns
    -------
    q_liq_top : float
        The specific flow rate of the liquid for 1 m of drain septum at top of column, [m**2/s]
    References
    ----------
    Дытнерский, страница 239, формула 6.39
    """
    return L_top / (rho_top_liq * L_septum)


@unitcheck(rho_bot_liq="kg/m**3", L_bot="kg/s", L_septum="m", res_unit="m**2/s")
def q_liq_bot(rho_bot_liq, L_septum, L_bot):
    """
    Calculates the specific flow rate of the liquid for 1 m of drain septum
    Parameters
    ----------
    rho_bot_liq : float
        The destiny of liquid at the bottom of column, [kg/m**3]
    L_bot : float
        The flow rate of liquid at the bottom of column, [kg/s]
    L_septum : float
        The length of drain septum, [m]
    Returns
    -------
    q_liq_top : float
        The specific flow rate of the liquid for 1 m of drain septum at the bottom of column, [m**2/s]
    References
    ----------
    Дытнерский, страница 239, формула 6.39
    """
    return L_bot / (rho_bot_liq * L_septum)


def epsi_vapor_top(Fr_top):
    """
    Calculates the vapor content of bubble layer at the top of column
    Parameters
    ----------
    Fr_top : float
        The Frudo criterion at the top of column, [dimensionless]
    Returns
    -------
    epsi_vapor_top : float
        The vapor content of bubble layer at the top of column, [dimensionless]
    References
    ----------
    Дытнерский, страница 207, формула 5.47
    """
    return Fr_top**0.5 / (1 + Fr_top**0.5)


def epsi_vapor_bot(Fr_bot):
    """
    Calculates the vapor content of bubble layer at the bottom of column
    Parameters
    ----------
    Fr_bot : float
        The Frudo criterion at the bottom of column, [dimensionless]
    Returns
    -------
    epsi_vapor_bot : float
        The vapor content of bubble layer at the bottom of column, [dimensionless]
    References
    ----------
    Дытнерский, страница 207, формула 5.47
    """
    return Fr_bot**0.5 / (1 + Fr_bot**0.5)  


@unitcheck(w_oper_top="m/s", heigth_layer_top="m", g="m/s**2")
def Fr_top(w_oper_top, heigth_layer_top):
    """
    Calculates the Frudo criterion at the top of column
    Parameters
    ----------
    w_oper_top : float
        The operating speed at the top of column, [m/s]
    heigth_layer_top : float
        The heigth ligth layer of  the liquid at the top of column, [m]
    Returns
    -------
    Fr_top : float
        The Frudo criterion at the top of column, [dimensionless]
    References
    ----------
    Дытнерский, страница 207, формула 5.47
    """
    return w_oper_top**2 / (g * heigth_layer_top)


@unitcheck(w_oper_bot="m/s", heigth_layer_bot="m", g="m/s**2")
def Fr_bot(w_oper_bot, heigth_layer_bot):
    """
    Calculates the Frudo criterion at the bottom of column
    Parameters
    ----------
    w_oper_bot: float
        The operating speed at the bottom of column, [m/s]
    heigth_layer_bot : float
        The heigth ligth layer of  the liquid at the bottom of column, [m]
    Returns
    -------
    Fr_bot : float
        The Frudo criterion at the bottom of column, [dimensionless]
    References
    ----------
    Дытнерский, страница 207, формула 5.47
    """
    return w_oper_bot**2 / (g * heigth_layer_bot)


@unitcheck(H_bwplate="m", h_bubble_top="m", res_unit="m")
def H_separate_top(H_bwplate, h_bubble_top):
    """
    Calculates the heigth of separation space at the top of column.
    Parameters
    ----------
    H_bwplate : float
        The heigth of between plates, [m]
    h_bubble : float
        The heigth of bubble layer at the top of column, [m]
    Returns
    -------
    H_separate_top : float
        The heigth of separation space at the top of column, [m]
    References
    ----------
    Дытнерский, страница 242, формула 6.42
    """    
    return H_bwplate - h_bubble_top


@unitcheck(H_bwplate="m", h_bubble_bot="m", res_unit="m")
def H_separate_bot(H_bwplate, h_bubble_bot):
    """
    Calculates the heigth of separation space at the bottom of column.
    Parameters
    ----------
    H_bwplate : float
        The heigth of between plates, [m]
    h_bubble : float
        The heigth of bubble layer at the bottom of column, [m]
    Returns
    -------
    H_separate_bot : float
        The heigth of separation space at the bottom of column, [m]
    References
    ----------
    Дытнерский, страница 242, формула 6.42
    """    
    return H_bwplate - h_bubble_bot


@unitcheck(heigth_layer_top="m", res_unit="m")
def h_bubble_top(heigth_layer_top, epsi_vapor_top):
    """
    Calculates the heigth of bubble layer at the top of column.
    Parameters
    ----------
    epsi_vapor_top : float
        The vapor content of bubble layer at the top of column, [dimensionless]
    heigth_layer_top : float
        The heigth ligth layer of  the liquid at the top of column, [m]
    Returns
    -------
    h_bubble_top : float
        The heigth of of bubble layer at the top of column, [m]
    References
    ----------
    Дытнерский, страница 242
    """
    return heigth_layer_top / (1 - epsi_vapor_top)


@unitcheck(heigth_layer_bot="m", res_unit="m")
def h_bubble_bot(heigth_layer_bot, epsi_vapor_bot):
    """
    Calculates the heigth of bubble layer at the bottom of column.
    Parameters
    ----------
    epsi_vapor_bot : float
        The vapor content of bubble layer at the bottom of column, [dimensionless]
    heigth_layer_bot : float
        The heigth ligth layer of  the liquid at the bottom of column, [m]
    Returns
    -------
    h_bubble_bot : float
        The heigth of of bubble layer at the bottom of column, [m]
    References
    ----------
    Дытнерский, страница 242
    """
    return heigth_layer_bot / (1 - epsi_vapor_bot)


@unitcheck(Massl="g/mol", Massh="g/mol", mu_solv="Pa*s", nul="sm**3/mol", nuh="sm**3/mol", res_unit="m**2/s")
def Diff_20(Massl, Massh, A , B, mu_solv, nul, nuh):
    """
    Calculates the diffusion coefficient at 20 degrees celcium.
    Parameters
    ----------
    Massl : float
        The molar mass of low-boilling component, [g/mol]
    Massh : float
        The molar mass of high-boilling component, [g/mol]
    A : float
        The correction coefficient depending on the properties solute, [dimensionless]
    B : float
        The correction coefficient depending on the properties solvent, [dimensionless]
    mu_solv : float
        The viscocity of solvent liquid, [Pa/s]
    nul : float
        The molar volume of solute, [sm**3/s]
    nuh : float
        The molar volume of solvent, [sm**3/s]
    Returns
    -------
    Diff_20 : float
        The diffusion coefficient at 20 degrees celcium, [m**2/s]
    References
    ----------
    Романков, страница 289, формула 6.22
    """
    return 1e-6 * ((1/Massl) + (1/Massh))**0.5 / (A * B * mu_solv**0.5 * ((nul)**0.66 + (nuh)*0.66)**2)


def b(mul_20, rhol_20):
    """
    Calculates the temperature coefficient.
    Parameters
    ----------
    mul_20 : float
        The viscocity of low-boilling component of liquid at 20 degrees celcium, [Pa/s]
    rhol_20 : float
        The destinity of low-boilling component of liquid at 20 degrees celcium, [kg/m**3]
    Returns
    -------
    b : float
        The temperature coefficient, [dimensionless]
    References
    ----------
    Романков, страница 289, формула 6.24
    """
    return mul_20**0.5/rhol_20**0.66


def Diff_liq(Diff_20, b, t_boil):
    """
    Calculates the diffusion coefficient of liquid phaze.
    Parameters
    ----------
    calc_Diffcoef20 : float
        The diffusion coefficient at 20 degrees celcium, [m**2/s]
    b : float
        The temperature coefficient, [dimensionless]
    t_boil : float
        The boiling temperature of liquid, [degrees celcium]
    Returns
    -------
    Diff_liq : float
        The diffusion coefficient of liquid phaze.
    References
    ----------
    Романков, страница 289, формула 6.23
    """
    return Diff_20 * (1 + b * (t_boil - 20))


@unitcheck(t_boil="K", Massl="g/mol", Massh="g/mol", P_abs="Pa", nul="sm**3/mol", nuh="sm**3/mol", res_unit="m**2/s")
def Diff_vapor(t_boil, P_abs, Massl, Massh, nul, nuh):
    """
    Calculates the diffusion coefficient of vapor.
    Parameters
    ----------
    Massl : float
        The molar mass of the low-boilling component, [g/mol]
    Massh : float
        The molar mass of the high-boilling component, [g/mol]
    t_boil : float
        The boiling temperature of the low-boiling component, [K]
    P : float
        The absolute pressure of the column, [Pa]
    nul : float
        The molar volume of solute, [sm**3/s]
    nuh : float
        The molar volume of solvent, [sm**3/s]
    Returns
    -------
    Diff_vapor : float
        The diffusion coefficient of vapor, [m**2/s]
    References
    ----------
    Романков, страница 234, формула 6.25
    """
    return 4.22e-2 * t_boil**(3/2) * ((1/Massl) + (1/Massh))**0.5 / (P_abs * ((nul)**0.66 + (nuh)*0.66)**2)


def beta_liq(Diff_liq, epsi_vapor, U_coef, heigth_layer, mu_vapor, mu_mix):
    """
    Calculates the coefficient masstransfer of liquid.
    Parameters
    ----------
    Diff_liq : float
        The diffusion coefficient of liquid phaze.
    epsi_vapor : float
        The vapor content of bubble layer, [dimensionless]
    heigth_layer : float
        The heigth ligth layer of  the liquid, [m]
    mu_mix : float
        The mix viscocity of liquid, [Pa/s]
    mu_vapor : float
        The mix viscocity of vapor, [Pa/s]
    U_coef : float
        The specific coefficient for beta_liq equation.
    Returns
    -------
    beta_liq : float
        The coefficient masstransfer of liquid, [m/s]
    References
    ----------
    Дытнерский, формула 6.37, стр.239
    """              
    return 6.24e+5 * (Diff_liq**0.5) * heigth_layer * ((U_coef/(1-epsi_vapor))**0.5) * (mu_vapor / (mu_vapor + mu_mix))**0.5


def beta_vapor(Diff_vapor, w_oper, epsi_vapor, heigth_layer, Fc, mu_vapor, mu_mix):
    """
    Calculates the coefficient masstransfer of vapor.
    Parameters
    ----------
    Diff_vapor : float
        The diffusion coefficient of vapor phaze, []
    epsi_vapor : float
        The vapor content of bubble layer, [dimensionless]
    heigth_layer : float
        The heigth ligth layer of  the liquid, [m]
    mu_mix : float
        The mix viscocity of liquid, [Pa/s]
    mu_vapor : float
        The mix viscocity of vapor, [Pa/s]
    Fc : float
        The free section of a plate, [dimensionless]
    Returns
    -------
    beta_vapor : float
        The coefficient masstransfer of vapor, [m/s]
    References
    ----------
    Дытнерский, формула 6.38, стр.239
    """
    return 6.24e+5 * Diff_vapor**0.5 * ((w_oper/epsi_vapor)**0.5) * heigth_layer * Fc * ((mu_vapor / (mu_vapor + mu_mix))**0.5)
#endregion