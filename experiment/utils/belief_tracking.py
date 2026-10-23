# failure_regions, force_osd_failures, low_weight_failures
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def pad_history(history: np.ndarray, max_iter: int) -> np.ndarray:
    history = np.asarray(history)
    n, m = history.shape
    if n >= max_iter:
        return history

    pad_rows = max_iter - n
    last_row = history[-1:]

    padding = np.repeat(last_row, pad_rows, axis=0)
    return np.vstack([history, padding])

# ---------------------------------------------------------------------------
# Belief Tracking Analysis
# ---------------------------------------------------------------------------

def analyze_belief_history(
    belief_history: np.ndarray,
    max_iter: int,
    pcm: np.ndarray,
    syndrome: np.ndarray,
    *,
    fixed_window: int = 3,
    fixed_delta_tol: float = 1e-9,
    oscillation_period: int = 2,
    oscillation_delta_tol: float = 1e-9,
) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    """
    Analyze BP belief / LLR history.

    Parameters
    ----------
    belief_history:
        Shape (num_iterations, num_variables).

    pcm:
        Parity-check matrix H, shape (num_checks, num_variables).

    syndrome:
        Target syndrome, shape (num_checks,).

    Returns
    -------
    summary_df:
        One-row global summary.

    variable_df:
        Per-variable summary, sorted by instability.

    info:
        Intermediate arrays for plotting/debugging.
    """
    belief_history = pad_history(belief_history, max_iter)
    pcm = np.asarray(pcm, dtype=np.uint8) % 2
    syndrome = np.asarray(syndrome, dtype=np.uint8) % 2

    num_iterations, num_variables = belief_history.shape

    if num_iterations < 2:
        raise ValueError("belief_history must contain at least 2 iterations.")

    if pcm.shape[1] != num_variables:
        raise ValueError(
            f"pcm has {pcm.shape[1]} variables, "
            f"but belief_history has {num_variables} variables."
        )

    if pcm.shape[0] != syndrome.shape[0]:
        raise ValueError(
            f"pcm has {pcm.shape[0]} checks, "
            f"but syndrome has length {syndrome.shape[0]}."
        )

    # ------------------------------------------------------------------
    # Basic histories
    # ------------------------------------------------------------------
    hard_history = (belief_history < 0).astype(np.uint8)

    belief_step = np.diff(belief_history, axis=0)
    hard_step = np.diff(hard_history, axis=0)

    llr_step_l2 = np.r_[np.nan, np.linalg.norm(belief_step, axis=1)]
    llr_step_linf = np.r_[np.nan, np.max(np.abs(belief_step), axis=1)]
    hard_step_hamming = np.r_[np.nan, np.count_nonzero(hard_step, axis=1)]

    syndrome_from_hard_history = (hard_history @ pcm.T) % 2
    syndrome_residual_history = syndrome_from_hard_history ^ syndrome[None, :]
    syndrome_residual_weight = syndrome_residual_history.sum(axis=1)

    sign_flip_count = np.count_nonzero(np.diff(hard_history, axis=0), axis=0)
    llr_range = np.ptp(belief_history, axis=0)
    final_abs_llr = np.abs(belief_history[-1])
    final_hard_decision = hard_history[-1]

    # ------------------------------------------------------------------
    # Fixed-point detection
    # ------------------------------------------------------------------
    if num_iterations >= fixed_window:
        tail_hard = hard_history[-fixed_window:]
        hard_is_fixed = np.all(tail_hard == tail_hard[-1])

        tail_step = np.diff(belief_history[-fixed_window:], axis=0)
        if tail_step.shape[0] == 0:
            tail_max_llr_step_l2 = 0.0
        else:
            tail_max_llr_step_l2 = float(
                np.max(np.linalg.norm(tail_step, axis=1))
            )

        fixed_point_detected = bool(
            hard_is_fixed and tail_max_llr_step_l2 <= fixed_delta_tol
        )
    else:
        tail_max_llr_step_l2 = np.nan
        fixed_point_detected = False

    # ------------------------------------------------------------------
    # Periodic oscillation detection
    # ------------------------------------------------------------------
    if num_iterations >= 2 * oscillation_period:
        previous_belief = belief_history[-2 * oscillation_period:-oscillation_period]
        current_belief = belief_history[-oscillation_period:]

        previous_hard = hard_history[-2 * oscillation_period:-oscillation_period]
        current_hard = hard_history[-oscillation_period:]

        oscillation_repeat_delta = float(np.max(np.linalg.norm(current_belief - previous_belief, axis=1)))

        repeats = np.array_equal(previous_hard, current_hard)
        nontrivial = any(
            not np.array_equal(current_hard[0], state)
            for state in current_hard[1:]
        )

        oscillation_detected = bool(
            repeats
            and nontrivial
            and oscillation_repeat_delta <= oscillation_delta_tol
        )
    else:
        oscillation_repeat_delta = np.nan
        oscillation_detected = False

    if fixed_point_detected:
        classification = f"Fixed Point(window={fixed_window})"
    elif oscillation_detected:
        classification = f"Oscillation(period={oscillation_period})"
    else:
        classification = "Unclassified"

    # ------------------------------------------------------------------
    # Summary tables
    # ------------------------------------------------------------------
    summary_df = pd.DataFrame(
        [
            {
                "num_iterations": num_iterations,
                "num_variables": num_variables,
                "num_detectors": pcm.shape[0],
                "classification": classification,
                "initial_syndrome_residual_weight": int(syndrome_residual_weight[0]),
                "final_syndrome_residual_weight": int(syndrome_residual_weight[-1]),
                "last_llr_step_l2": float(llr_step_l2[-1]) if not np.isnan(llr_step_l2[-1]) else np.nan,
                "last_llr_step_linf": float(llr_step_linf[-1]) if not np.isnan(llr_step_linf[-1]) else np.nan,
                "tail_max_llr_step_l2": tail_max_llr_step_l2,
                "oscillation_repeat_delta": oscillation_repeat_delta,
                "last_hard_step_hamming": float(hard_step_hamming[-1]) if not np.isnan(hard_step_hamming[-1]) else np.nan,
                "total_sign_flips": int(sign_flip_count.sum()),
            }
        ]
    )

    variable_df = pd.DataFrame(
        {
            "variable": np.arange(num_variables),
            "llr_initial": belief_history[0],
            "llr_final": belief_history[-1],
            "abs_llr_final": final_abs_llr,
            "llr_range": llr_range,
            "sign_flip_count": sign_flip_count,
            "final_hard_decision": final_hard_decision,
        }
    )

    info = {
        "belief_history": belief_history,
        "hard_history": hard_history,
        "belief_step": belief_step,
        "hard_step": hard_step,
        "llr_step_l2": llr_step_l2,
        "llr_step_linf": llr_step_linf,
        "hard_step_hamming": hard_step_hamming,
        "syndrome_from_hard_history": syndrome_from_hard_history,
        "syndrome_residual_history": syndrome_residual_history,
        "syndrome_residual_weight": syndrome_residual_weight,
        "classification": classification,
        "tail_max_llr_step_l2": tail_max_llr_step_l2,
        "oscillation_repeat_delta": oscillation_repeat_delta,
    }

    return summary_df, variable_df, info


# ---------------------------------------------------------------------------
# Plot Belief Tracking
# ---------------------------------------------------------------------------

def plot_belief_tracking(
    belief_history: np.ndarray,
    ax: plt.Axes
):
    """
    Plot the belief history of BP. 
    
    Parameters:
    -------------
    belief_history: np.ndarray
        A 2D array of shape (num_iterations, num_variables), where each entry is the belief (LLR) of a variable at a given iteration.
    """
    belief_delta = np.ptp(belief_history, axis=0)
    variables_to_plot = np.argsort(belief_delta)[::-1]

    for var in variables_to_plot:
        ax.plot(belief_history[:, var], marker='o', label=f'var {var}')

    ax.axhline(0, color='black', linewidth=1, alpha=0.4)
    ax.set_xlabel('iteration')
    ax.set_ylabel('belief / LLR')
    ax.legend(ncol=2, fontsize=8)
    ax.grid(True, alpha=0.3)


def plot_history_heatmap(
    history: np.ndarray,
    kind: str,
    fig: plt.Figure,
    ax: plt.Axes
):
    """
    Plot the belief history of BP as a heatmap. 
    
    Parameters:
    -------------
    history: np.ndarray
        A 2D array of shape (num_iterations, num_variables).
    kind: str
        The type of data to plot (e.g., 'belief', 'hard decision').
    """
    im = ax.imshow(history.T, aspect='auto', origin='lower')
    ax.set_xlabel('iteration')
    ax.set_ylabel('variable index')
    fig.colorbar(im, ax=ax, label=kind)


def plot_belief_analysis(
    info: dict,
    variable_df: pd.DataFrame,
    *,
    top_k: int = 8,
    variables: list[int] | None = None,
    figsize: tuple = (14, 8),
):
    """
    Plot belief dynamics.

    If `variables` is None, plot the top-k most unstable variables.
    Otherwise, plot the specified variables.
    """
    belief_history = info["belief_history"]
    llr_step_l2 = info["llr_step_l2"]
    llr_step_linf = info["llr_step_linf"]
    hard_step_hamming = info["hard_step_hamming"]
    syndrome_residual_weight = info["syndrome_residual_weight"]
    classification = info["classification"]

    num_iterations = belief_history.shape[0]
    iteration_axis = np.arange(num_iterations)
    
    if variables is None:
        variable_df = variable_df.sort_values(
            ["sign_flip_count", "llr_range", "abs_llr_final"],
            ascending=[False, False, True],
        ).reset_index(drop=True)
        variables = variable_df.head(top_k)["variable"].to_numpy(dtype=int)
    else:
        variables = np.asarray(variables, dtype=int)

    fig, axes = plt.subplots(2, 2, figsize=figsize)

    axes[0, 0].plot(
        iteration_axis,
        np.nan_to_num(llr_step_l2, nan=0.0),
        marker="o",
        label="LLR step L2",
    )
    axes[0, 0].plot(
        iteration_axis,
        np.nan_to_num(llr_step_linf, nan=0.0),
        marker="s",
        label="LLR step Linf",
    )
    axes[0, 0].set_xlabel("iteration")
    axes[0, 0].set_ylabel("step norm")
    axes[0, 0].set_title("LLR update size")
    axes[0, 0].legend()

    axes[0, 1].plot(
        iteration_axis,
        syndrome_residual_weight,
        marker="o",
        label="syndrome residual weight",
    )
    axes[0, 1].plot(
        iteration_axis,
        np.nan_to_num(hard_step_hamming, nan=0.0),
        marker="s",
        label="hard-decision Hamming step",
    )
    axes[0, 1].set_xlabel("iteration")
    axes[0, 1].set_ylabel("count")
    axes[0, 1].set_title("Hard decision and syndrome residual")
    axes[0, 1].legend()

    for variable in variables:
        axes[1, 0].plot(
            iteration_axis,
            belief_history[:, variable],
            marker="o",
            label=f"var {variable}",
        )

    axes[1, 0].axhline(0, color="black", linewidth=1, linestyle="--")
    axes[1, 0].set_xlabel("iteration")
    axes[1, 0].set_ylabel("belief / LLR")
    axes[1, 0].set_title("Selected variable trajectories")
    axes[1, 0].legend(ncol=2, fontsize=8)

    selected_df = variable_df.set_index("variable").loc[variables].reset_index()

    axes[1, 1].bar(
        np.arange(len(variables)),
        selected_df["sign_flip_count"],
    )
    axes[1, 1].set_xticks(np.arange(len(variables)))
    axes[1, 1].set_xticklabels([str(v) for v in variables])
    axes[1, 1].set_xlabel("variable")
    axes[1, 1].set_ylabel("sign flip count")
    axes[1, 1].set_title("Variable sign instability")

    fig.suptitle(f"Belief dynamics: {classification}", fontsize=14)
    fig.tight_layout()

    return fig, axes

