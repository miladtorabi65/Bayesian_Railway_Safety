from manim import *


class RailwaySafetyModels(MovingCameraScene):
    """Animated explanation of the progression M1 -> M2 -> M3.

    Render a quick preview:
        manim -pql railway_safety_models_manim.py RailwaySafetyModels

    Render in high quality:
        manim -pqh railway_safety_models_manim.py RailwaySafetyModels
    """

    BG = "#0E1117"
    PANEL = "#171C26"

    def construct(self):
        self.camera.background_color = self.BG
        self.intro()
        self.notation()
        self.common_structure()
        self.log_link()
        self.offset()
        self.m1_poisson()
        self.m2_negative_binomial()
        self.m3_hierarchical()
        self.hierarchy_and_effects()
        self.sum_to_zero()
        self.summary()

    # ---------- helpers ----------
    def clear(self):
        self.play(*[FadeOut(m) for m in self.mobjects], run_time=0.7)

    def heading(self, text):
        h = Text(text, font_size=43, weight=BOLD, color=BLUE_B)
        h.to_edge(UP)
        return h

    def caption(self, text):
        c = Text(text, font_size=26, color=GREY_A, line_spacing=0.9)
        c.to_edge(DOWN, buff=0.35)
        return c

    def section(self, text):
        box = RoundedRectangle(
            width=11.4, height=2.1, corner_radius=0.22,
            fill_color=self.PANEL, fill_opacity=1,
            stroke_color=BLUE_D, stroke_width=2,
        )
        label = Text(text, font_size=46, weight=BOLD)
        label.move_to(box)
        self.play(FadeIn(box), Write(label))
        self.wait(1.0)
        self.play(FadeOut(VGroup(box, label)))

    def info_box(self, symbol, text, color):
        s = MathTex(symbol, font_size=43, color=color)
        t = Text(text, font_size=24, color=WHITE)
        row = VGroup(s, t).arrange(RIGHT, buff=0.3)
        frame = SurroundingRectangle(row, color=color, buff=0.18, corner_radius=0.12)
        return VGroup(frame, row)

    # ---------- scenes ----------
    def intro(self):
        title = Text(
            "From Poisson to a Hierarchical\nNegative Binomial Model",
            font_size=51, weight=BOLD, line_spacing=0.9,
        )
        sub = Text(
            "Significant railway accidents across European countries",
            font_size=29, color=BLUE_B,
        ).next_to(title, DOWN, buff=0.45)
        path = MathTex(
            r"M_1:\ \text{Poisson}", r"\longrightarrow",
            r"M_2:\ \text{Negative Binomial}", r"\longrightarrow",
            r"M_3:\ \text{Hierarchical NB}", font_size=34,
        ).next_to(sub, DOWN, buff=0.65)
        path[0].set_color(BLUE_C)
        path[2].set_color(YELLOW_C)
        path[4].set_color(GREEN_C)
        self.play(Write(title), run_time=1.4)
        self.play(FadeIn(sub, shift=UP * 0.2))
        self.play(Write(path), run_time=1.4)
        self.wait(2)
        self.clear()

    def notation(self):
        self.section("Data and notation")
        h = self.heading("One observation is one country in one year")
        self.play(Write(h))

        table = MathTable(
            [["Italy", "2010", "Y", "E"],
             ["Italy", "2011", "Y", "E"],
             [r"\vdots", r"\vdots", r"\vdots", r"\vdots"],
             ["France", "2024", "Y", "E"]],
            col_labels=[Text("Country", font_size=22), Text("Year", font_size=22),
                        Text("Accidents", font_size=22), Text("Train-km", font_size=22)],
            include_outer_lines=True,
            line_config={"stroke_color": GREY_C},
        ).scale(0.67).shift(LEFT * 3.1 + DOWN * 0.2)

        boxes = VGroup(
            self.info_box(r"Y_{it}", "observed accident count", RED_C),
            self.info_box(r"E_{it}", "train-km exposure", BLUE_C),
            self.info_box(r"x_t", "centred year = year - 2017", YELLOW_C),
            self.info_box(r"\mu_{it}", "expected accident count", TEAL_C),
        ).arrange(DOWN, aligned_edge=LEFT, buff=0.28).scale(0.85)
        boxes.shift(RIGHT * 3.05 + DOWN * 0.15)

        idx = MathTex(r"i=\text{country},\qquad t=\text{year}", font_size=32, color=GREY_A)
        idx.next_to(boxes, DOWN, buff=0.35)

        self.play(Create(table), run_time=1.2)
        self.play(LaggedStart(*[FadeIn(b, shift=LEFT * 0.15) for b in boxes], lag_ratio=0.2))
        self.play(Write(idx))
        self.play(FadeIn(self.caption("All three models use the same 405 country-year observations.")))
        self.wait(2.3)
        self.clear()

    def common_structure(self):
        self.section("The common mean structure")
        h = self.heading("Expected accidents = exposure × underlying rate")
        self.play(Write(h))
        eq = MathTex(
            r"\log(\mu_{it})", "=", r"\log(E_{it})", "+", r"\alpha", "+", r"\beta x_t",
            font_size=57,
        ).shift(UP * 0.8)
        eq[0].set_color(RED_B)
        eq[2].set_color(BLUE_C)
        eq[4].set_color(PURPLE_B)
        eq[6].set_color(YELLOW_C)
        self.play(Write(eq))

        labels = ["expected count", "exposure offset", "baseline log-rate", "time trend"]
        parts = [eq[0], eq[2], eq[4], eq[6]]
        colors = [RED_B, BLUE_C, PURPLE_B, YELLOW_C]
        annotations = VGroup()
        for part, label, color in zip(parts, labels, colors):
            brace = Brace(part, DOWN, color=color)
            txt = Text(label, font_size=20, color=color).next_to(brace, DOWN, buff=0.1)
            annotations.add(VGroup(brace, txt))
        self.play(LaggedStart(*[FadeIn(a) for a in annotations], lag_ratio=0.2))

        exp_eq = MathTex(
            r"\mu_{it}=E_{it}\times\exp(\alpha+\beta x_t)", font_size=52,
        ).shift(DOWN * 1.65)
        exp_eq.set_color_by_tex(r"\mu_{it}", RED_B)
        exp_eq.set_color_by_tex(r"E_{it}", BLUE_C)
        self.play(TransformFromCopy(eq, exp_eq))
        self.play(FadeIn(self.caption("The model predicts a count by multiplying exposure by an underlying accident rate.")))
        self.wait(2.2)
        self.clear()

    def log_link(self):
        self.section("Why use a logarithm?")
        h = self.heading("The expected accident count must always be positive")
        self.play(Write(h))
        bad = MathTex(r"\mu_{it}=\alpha+\beta x_t", font_size=52, color=RED_C).shift(UP * 1)
        bad_note = MathTex(r"\text{could give }\mu_{it}<0", font_size=39, color=RED_B).next_to(bad, DOWN)
        self.play(Write(bad), Write(bad_note))
        self.play(Create(Cross(VGroup(bad, bad_note), stroke_color=RED_C, stroke_width=7)))
        good = MathTex(
            r"\log(\mu_{it})=\eta_{it}", r"\Longrightarrow", r"\mu_{it}=\exp(\eta_{it})>0",
            font_size=50,
        ).shift(DOWN * 1.25)
        good[0].set_color(BLUE_C)
        good[2].set_color(GREEN_C)
        self.play(Write(good))
        self.play(FadeIn(self.caption("Exponentiation guarantees a positive expected count.")))
        self.wait(2)
        self.clear()

    def offset(self):
        self.section("Why use train-km as an offset?")
        h = self.heading("More traffic means more opportunities for accidents")
        self.play(Write(h))

        a = VGroup(Text("Country A", font_size=30, weight=BOLD),
                   MathTex(r"E_A=100", font_size=39, color=BLUE_C),
                   MathTex(r"\mu_A=10", font_size=39, color=RED_C)).arrange(DOWN, buff=0.25)
        b = VGroup(Text("Country B", font_size=30, weight=BOLD),
                   MathTex(r"E_B=200", font_size=39, color=BLUE_C),
                   MathTex(r"\mu_B=20", font_size=39, color=RED_C)).arrange(DOWN, buff=0.25)
        ga = VGroup(SurroundingRectangle(a, color=BLUE_D, buff=0.3), a).shift(LEFT * 3.2 + UP * 0.25)
        gb = VGroup(SurroundingRectangle(b, color=BLUE_D, buff=0.3), b).shift(RIGHT * 3.2 + UP * 0.25)
        ar = Arrow(ga.get_right(), gb.get_left(), color=YELLOW_C, buff=0.3)
        note = Text("twice the exposure\n→ twice the expected count", font_size=26, color=YELLOW_C)
        note.next_to(ar, UP, buff=0.15)
        eq = MathTex(r"\mu_{it}=E_{it}\times\exp(\alpha+\beta x_t)", font_size=50).shift(DOWN * 1.8)
        fixed = Text("The offset coefficient is fixed at 1; it is not estimated.", font_size=26, color=GREY_A)
        fixed.next_to(eq, DOWN, buff=0.35)
        self.play(FadeIn(ga), FadeIn(gb))
        self.play(GrowArrow(ar), FadeIn(note))
        self.play(Write(eq), FadeIn(fixed))
        self.wait(2.2)
        self.clear()

    def m1_poisson(self):
        self.section("Model 1 — Poisson regression")
        h = self.heading("M1: Exposure + time trend + ordinary count noise")
        self.play(Write(h))
        like = MathTex(r"Y_{it}\sim\operatorname{Poisson}(\mu_{it})", font_size=56).shift(UP * 1.45)
        mean = MathTex(r"\log(\mu_{it})=\log(E_{it})+\alpha+\beta x_t", font_size=50).shift(UP * 0.35)
        mean.set_color_by_tex(r"E_{it}", BLUE_C)
        mean.set_color_by_tex(r"\alpha", PURPLE_B)
        mean.set_color_by_tex(r"\beta", YELLOW_C)
        var = MathTex(r"\operatorname{Var}(Y_{it}\mid\mu_{it})=\mu_{it}", font_size=47, color=TEAL_B).shift(DOWN * 0.9)
        q = Text("Question: Is ordinary Poisson randomness enough?", font_size=30)
        qbox = SurroundingRectangle(q, color=BLUE_C, buff=0.28, corner_radius=0.14)
        qg = VGroup(qbox, q).shift(DOWN * 2.15)
        self.play(Write(like), Write(mean), Write(var))
        self.play(FadeIn(qg))
        self.wait(1.8)

        self.play(FadeOut(like), FadeOut(var), FadeOut(qg), mean.animate.shift(UP * 0.9))
        alpha = self.info_box(r"\alpha", "baseline log accident rate in 2017", PURPLE_B).shift(LEFT * 2.5 + DOWN * 0.25)
        center = MathTex(r"x_t=\text{year}-2017,\qquad x_{2017}=0", font_size=39, color=YELLOW_C).shift(RIGHT * 2.2 + DOWN * 0.25)
        self.play(FadeIn(alpha), Write(center))
        baseline = MathTex(r"\exp(\alpha)=\text{baseline accident rate}", font_size=42).shift(DOWN * 1.65)
        self.play(Write(baseline))
        self.wait(1.6)
        self.play(FadeOut(alpha), FadeOut(center), FadeOut(baseline))

        beta = self.info_box(r"\beta", "common yearly change on the log-rate scale", YELLOW_C).shift(UP * 0.1)
        rr = MathTex(r"\exp(\beta)=\text{annual rate ratio}", font_size=44).next_to(beta, DOWN, buff=0.45)
        ex = MathTex(r"\exp(\beta)=0.97\Longrightarrow100(0.97-1)=-3\%", font_size=41).next_to(rr, DOWN, buff=0.5)
        ex.set_color_by_tex("0.97", GREEN_C)
        ex.set_color_by_tex("-3", GREEN_C)
        sentence = Text("The accident rate decreases by 3% per year.", font_size=28, color=GREEN_C).next_to(ex, DOWN, buff=0.35)
        self.play(FadeIn(beta), Write(rr), Write(ex), FadeIn(sentence))
        self.wait(2.2)
        self.clear()

    def m2_negative_binomial(self):
        self.section("Model 2 — Negative binomial regression")
        h = self.heading("M2: Keep the same mean, but allow extra variation")
        self.play(Write(h))

        mean = MathTex(r"\log(\mu_{it})=\log(E_{it})+\alpha+\beta x_t", font_size=47).shift(UP * 1.55)
        unchanged = Text("unchanged mean structure", font_size=24, color=GREEN_C).next_to(mean, UP, buff=0.15)
        like = MathTex(r"Y_{it}\sim\operatorname{NB}(\mu_{it},r)", font_size=53).shift(UP * 0.35)
        var = MathTex(r"\operatorname{Var}(Y_{it})=\mu_{it}+\frac{\mu_{it}^2}{r}", font_size=54).shift(DOWN * 0.9)
        var.set_color_by_tex(r"\frac{\mu_{it}^2}{r}", YELLOW_C)
        self.play(FadeIn(unchanged), Write(mean), Write(like), Write(var))
        brace = Brace(var[-1], DOWN, color=YELLOW_C)
        label = Text("extra-Poisson variation", font_size=24, color=YELLOW_C).next_to(brace, DOWN, buff=0.1)
        self.play(GrowFromCenter(brace), FadeIn(label))
        self.wait(1.6)

        self.play(FadeOut(unchanged), FadeOut(mean), FadeOut(like), FadeOut(brace), FadeOut(label), var.animate.shift(UP * 1.7))
        toy = Text("Toy example: both models expect 50 accidents", font_size=30).next_to(var, DOWN, buff=0.55)
        p = MathTex(r"\text{Poisson: }\operatorname{Var}(Y)=50", font_size=41, color=BLUE_C).next_to(toy, DOWN, buff=0.4)
        nb = MathTex(r"\text{NB with }r=10:\quad 50+\frac{50^2}{10}=300", font_size=41, color=YELLOW_C).next_to(p, DOWN, buff=0.35)
        self.play(FadeIn(toy), Write(p), Write(nb))
        self.wait(1.8)
        self.play(FadeOut(toy), FadeOut(p), FadeOut(nb))

        rbox = self.info_box(r"r", "negative-binomial dispersion parameter", YELLOW_C).shift(UP * 0.05)
        small = MathTex(r"\text{small }r\Rightarrow\text{strong overdispersion}", font_size=37, color=RED_B).next_to(rbox, DOWN, buff=0.45)
        large = MathTex(r"\text{large }r\Rightarrow\operatorname{Var}(Y)\approx\mu", font_size=37, color=GREEN_C).next_to(small, DOWN, buff=0.35)
        self.play(FadeIn(rbox), Write(small), Write(large))
        self.play(FadeIn(self.caption("M1 → M2 asks whether extra random variation is needed.")))
        self.wait(2.2)
        self.clear()

    def m3_hierarchical(self):
        self.section("Model 3 — Hierarchical negative binomial")
        h = self.heading("M3: Add persistent country-specific differences")
        self.play(Write(h))

        like = MathTex(r"Y_{it}\sim\operatorname{NB}(\mu_{it},r)", font_size=49).shift(UP * 1.5)
        base = MathTex(r"\log(\mu_{it})=\log(E_{it})+\alpha+\beta x_t", font_size=48).shift(UP * 0.35 + LEFT * 0.6)
        ui = MathTex(r"+u_i", font_size=50, color=GREEN_C).next_to(base, RIGHT, buff=0.15)
        prior = MathTex(r"u_i\sim\mathcal N(0,\sigma_{\mathrm{country}}^2)", font_size=48).shift(DOWN * 0.95)
        prior.set_color_by_tex(r"u_i", GREEN_C)
        prior.set_color_by_tex(r"\sigma_{\mathrm{country}}", PURPLE_B)
        self.play(Write(like), Write(base))
        self.play(GrowFromCenter(ui))
        self.play(ui.animate.scale(1.12), run_time=0.25)
        self.play(ui.animate.scale(1 / 1.12), run_time=0.25)
        self.play(Write(prior))
        q = Text("Do countries still differ after exposure, time, and overdispersion?", font_size=28)
        qg = VGroup(SurroundingRectangle(q, color=GREEN_C, buff=0.25), q).shift(DOWN * 2.15)
        self.play(FadeIn(qg))
        self.wait(2.2)
        self.clear()

    def hierarchy_and_effects(self):
        self.section("How the hierarchy works")
        h = self.heading("All 15 years from one country share the same effect")
        self.play(Write(h))

        centers = [LEFT * 4 + UP * 1.3, UP * 1.3, RIGHT * 4 + UP * 1.3]
        names = ["Country 1", "Country 2", "Country 27"]
        colors = [BLUE_C, YELLOW_C, GREEN_C]
        groups = VGroup()
        effects = VGroup()
        for center, name, color, idx in zip(centers, names, colors, ["1", "2", "27"]):
            circle = Circle(radius=0.82, color=color, fill_color=color, fill_opacity=0.12).move_to(center)
            label = Text(name, font_size=24, color=color).move_to(circle)
            effect = MathTex(rf"u_{{{idx}}}", font_size=38, color=color).next_to(circle, UP, buff=0.22)
            years = VGroup()
            lines = VGroup()
            for j, year in enumerate(["2010", "2011", "⋯", "2024"]):
                dot = Circle(radius=0.22, color=color, fill_color=color, fill_opacity=0.65)
                txt = Text(year, font_size=14).move_to(dot)
                item = VGroup(dot, txt).move_to(center + DOWN * 2 + RIGHT * (j - 1.5) * 0.64)
                years.add(item)
                lines.add(Line(circle.get_bottom(), item.get_top(), color=color, stroke_opacity=0.55))
            groups.add(VGroup(lines, circle, label, years))
            effects.add(effect)
        self.play(LaggedStart(*[FadeIn(g) for g in groups], lag_ratio=0.25), run_time=1.5)
        self.play(LaggedStart(*[Write(e) for e in effects], lag_ratio=0.25))
        self.play(FadeIn(self.caption("M2 models extra noise; M3 also learns persistent country structure.")))
        self.wait(2.1)
        self.clear()

        self.section("Interpreting country effects")
        h = self.heading("The relative country effect is exp(uᵢ)")
        self.play(Write(h))
        pos = VGroup(MathTex(r"u_i=0.4", font_size=45, color=GREEN_C),
                     MathTex(r"\exp(0.4)\approx1.49", font_size=43, color=GREEN_B),
                     Text("about 49% above the reference rate", font_size=25)).arrange(DOWN, buff=0.3).shift(LEFT * 3 + DOWN * 0.2)
        neg = VGroup(MathTex(r"u_i=-0.4", font_size=45, color=BLUE_C),
                     MathTex(r"\exp(-0.4)\approx0.67", font_size=43, color=BLUE_B),
                     Text("about 33% below the reference rate", font_size=25)).arrange(DOWN, buff=0.3).shift(RIGHT * 3 + DOWN * 0.2)
        self.play(FadeIn(pos), Create(Line(UP * 1.4, DOWN * 1.8, color=GREY_C)), FadeIn(neg))
        sigma = self.info_box(r"\sigma_{\mathrm{country}}", "spread of persistent country effects", PURPLE_B)
        sigma.to_edge(DOWN, buff=0.35)
        self.play(FadeIn(sigma))
        self.wait(2.2)
        self.clear()

    def sum_to_zero(self):
        self.section("Why use sum-to-zero country effects?")
        h = self.heading("Separating the overall intercept from the country effects")
        self.play(Write(h))
        a = MathTex(r"\alpha+u_i=-0.7+0.4=-0.3", font_size=46).shift(UP * 1.15)
        b = MathTex(r"\alpha+u_i=-0.9+0.6=-0.3", font_size=46).shift(UP * 0.2)
        brace = Brace(VGroup(a, b), RIGHT, color=RED_C)
        same = Text("same prediction", font_size=25, color=RED_C).next_to(brace, RIGHT, buff=0.15)
        problem = Text(
            "The intercept can move down while every country effect moves up.\n"
            "This creates redundant movement and slower MCMC mixing.",
            font_size=27, color=GREY_A, line_spacing=0.9,
        ).shift(DOWN * 1.1)
        self.play(Write(a), Write(b), GrowFromCenter(brace), FadeIn(same), FadeIn(problem))
        self.wait(1.8)
        self.play(FadeOut(VGroup(a, b, brace, same, problem)))

        constraint = MathTex(r"\boxed{\sum_{i=1}^{27}u_i=0}", font_size=61, color=GREEN_C).shift(UP * 1.0)
        example = MathTex(r"u_1=0.5,\qquad u_2=-0.2,\qquad u_3=-0.3", font_size=41).shift(DOWN * 0.1)
        total = MathTex(r"0.5-0.2-0.3=0", font_size=43, color=GREEN_C).shift(DOWN * 1.05)
        benefits = VGroup(
            Text("α = common overall reference level", font_size=26),
            Text("uᵢ = country deviations around that reference", font_size=26),
            Text("MCMC sampling becomes more stable", font_size=26),
        ).arrange(DOWN, aligned_edge=LEFT, buff=0.16).to_edge(DOWN, buff=0.28)
        self.play(Write(constraint), Write(example), Write(total), FadeIn(benefits))
        self.wait(2.5)
        self.clear()

    def summary(self):
        self.section("The complete model progression")
        h = self.heading("Each model solves a limitation of the previous one")
        self.play(Write(h))

        frames = []
        specs = [
            (LEFT * 4.1, BLUE_C, "M1", "Poisson", r"\operatorname{Var}(Y)=\mu",
             "Exposure + time trend\n+ ordinary count noise", "Is Poisson variation enough?"),
            (ORIGIN, YELLOW_C, "M2", "Negative binomial", r"\operatorname{Var}(Y)=\mu+\mu^2/r",
             "Same mean structure\n+ extra random variation", "Is overdispersion needed?"),
            (RIGHT * 4.1, GREEN_C, "M3", "Hierarchical NB", r"+u_i",
             "Extra variation\n+ persistent country effects", "Do countries still differ?"),
        ]
        for pos, color, code, name, formula, body, question in specs:
            rect = RoundedRectangle(width=3.55, height=4.25, corner_radius=0.2,
                                    fill_color=self.PANEL, fill_opacity=1,
                                    stroke_color=color, stroke_width=3).move_to(pos + DOWN * 0.15)
            content = VGroup(
                Text(code, font_size=36, weight=BOLD, color=color),
                Text(name, font_size=27, weight=BOLD),
                MathTex(formula, font_size=29, color=color),
                Text(body, font_size=22, line_spacing=0.9),
                Text(question, font_size=20, color=color),
            ).arrange(DOWN, buff=0.28).move_to(rect)
            frames.append(VGroup(rect, content))

        ar1 = Arrow(frames[0].get_right(), frames[1].get_left(), buff=0.15, color=GREY_A)
        ar2 = Arrow(frames[1].get_right(), frames[2].get_left(), buff=0.15, color=GREY_A)
        self.play(FadeIn(frames[0]))
        self.play(GrowArrow(ar1), FadeIn(frames[1]))
        self.play(GrowArrow(ar2), FadeIn(frames[2]))
        final = MathTex(
            r"\text{ordinary variation}\longrightarrow"
            r"\text{extra variation}\longrightarrow"
            r"\text{persistent country heterogeneity}",
            font_size=32,
        ).to_edge(DOWN, buff=0.25)
        final.set_color_by_tex("ordinary", BLUE_C)
        final.set_color_by_tex("extra", YELLOW_C)
        final.set_color_by_tex("persistent", GREEN_C)
        self.play(Write(final))
        self.wait(3)
