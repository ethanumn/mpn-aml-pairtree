import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

F_SIZE = (11.5,8)

plt.rcParams["figure.figsize"] = F_SIZE


class MPN_AML_METRICS_PDF():

    def __init__(self, filename="example.metrics.pdf", title="", details="", toc=True):

        self.file_name = filename
        self.title = title
        self.details = details

        self.figs = []
        self.toc = toc
        self.toc_text = []


    def render(self):

        with PdfPages(self.file_name) as pp:

            self._render_title_page(pp, self.title, self.details)
            self._render_toc(pp)

            for fig in self.figs:
                pp.savefig(fig)
                plt.close()


    def _render_title_page(self, pp, title, details):

        title_page = plt.figure()
        plt.axis("off")

        plt.text(0.5, 0.9, "%s" % title, ha="center", va="top", fontsize=32)
        plt.text(0.5, 0.8, "%s" % details, ha="center", va="top", fontsize=12)

        pp.savefig(title_page)
        plt.close()


    def _render_toc(self, pp):

        toc = plt.figure()
        plt.axis("off")
        plt.text(0.0, 0.8, "\n".join(self.toc_text), va="center", fontsize=12, transform = toc.axes[0].transAxes)

        pp.savefig(toc, bbox_inches='tight')
        plt.close()


    def add_details(self, conditions=[], statements=[], title=""):

        plot = plt.figure()
        plt.axis("off")
        y_pos = 1.0

        plt.text(0.5, y_pos, title, ha="center", va="center", fontsize=16)
        y_pos -= 0.1

        if len(conditions) == len(statements):

            for cond, text in zip(conditions, statements):
                y_pos -= 0.05
                plt.text(0.5, y_pos, ("PASSED: %s\n" if cond else "FAILED: %s\n") % text, ha="center", va="center", fontsize=8)

        self.add_figure(plot, title)
        plt.close()




    def add_figure(self, figure, title=""):

        plt.text(1.1, 1.1, "%d" % (len(self.figs) + 1), fontsize=12, transform = figure.axes[0].transAxes)

        if title:
            self.toc_text.append("%-8d%s" % (len(self.toc_text) + 2, title))
        else:
            self.toc_text.append("")

        self.figs.append(figure)
        plt.close()



    def add_table(self, data, title="",
                  header_color='#40466e', row_colors=['#f1f1f2', 'w'], edge_color='w',
                  bbox=[0, 0, 1, 1], header_columns=0, tight_layout=True):

        """
        Render a matplotlib table and add it to list of figures
        """

        fig, ax = plt.subplots()
        ax.axis('off')

        mpl_table = ax.table(cellText=data.values, bbox=bbox, colLabels=data.columns, loc="center")
        mpl_table.auto_set_font_size(True)

        for k, cell in mpl_table._cells.items():
            cell.set_edgecolor(edge_color)
            if k[0] == 0 or k[1] < header_columns:
                cell.set_text_props(weight='bold', color='w')
                cell.set_facecolor(header_color)
            else:
                cell.set_facecolor(row_colors[k[0]%len(row_colors) ])

        plt.suptitle(title)

        if tight_layout:
            plt.tight_layout()

        self.add_figure(fig, title)
        plt.close()


    def add_plot(self, x_data, y_data,
                 suptitle="", title="",
                 xlabel="", ylabel="",
                 xtick_rot=0, ylim=None,
                 color=["blue", "red"],
                 caption=""):
        """
        Render a matplotlib bar graph and add it to list of figures
        """

        plot = plt.figure()

        plt.bar(x_data, y_data, color=color)
        plt.xlabel(xlabel)
        plt.xticks(rotation=xtick_rot, fontsize=plot.axes[0].get_xticklabels()[0].get_fontsize() / 2) # workaround to make text fit
        plt.ylabel(ylabel)
        plt.suptitle(suptitle)
        plt.title(title)

        if caption:
            plt.text(0.15, 0.94, caption, ha='center', va='center', transform = plot.axes[0].transAxes)

        plt.tight_layout()


        if ylim:
            plt.ylim(ylim)

        self.add_figure(plot, "%s, %s" % (title, suptitle))
        plt.close()
