import Pkg
Pkg.add("PyCall")
Pkg.add("PyPlot")
using PyCall
np = pyimport("numpy")
using PyPlot: plt


mutable struct Solver_
    nu      # параметр, имеющий физический смысл в уравнении
    L       # длина отрезка интегрирования по х
    N       # количество узлов расчетной сетки
    t0      # начальный момент времени
    tN      # конечный момент времени
    dt      # шаг по времени
    k       # коэффициенты k для вычисления операторов дифференцирования
    FL      # оператор дифференцирования для линейной части
    FN      # оператор дифференцирования для нелинейной части
    t       # сетка по времени
    x       # сетка по пространству
    u       # решение
    solver  # метод класса, в котором реализован решатель
    show    # метод класса, демонстрирующий решение

    function Solver_(nu, L, N, t0, tN, dt)
        this = new()

        this.nu = nu
        this.L = L
        this.N = N
        this.t0 = t0
        this.tN = tN
        this.dt = dt

        this.k = np.arange(-N/2, N/2, 1)
        this.FL = (((2 * np.pi) / L ) * this.k) .^ 2 - this.nu * (((2 * np.pi) / L ) * this.k) .^ 4
        this.FN = - (1 / 2) * ((1im) * ((2 * np.pi) / L ) * this.k)

        this.solver = function()
            # количество шагов по времени
            nt = Int((this.tN - this.t0) /this.dt) 

            # сетка
            this.t = np.linspace(start=this.t0, stop=this.tN, num=nt)
            this.x = np.linspace(start=0, stop=this.L, num=this.N) 

            # начальные условия
            u0 = -np.cos((2 * np.pi * this.x) / this.L) + np.sin((4 * np.pi * this.x) / this.L)
            # коэффициенты пространства Фурье
            u0_hat = (1 / this.N) * np.fft.fftshift(np.fft.fft(u0))
            # отдельно для нелинейной части
            u0_hat2 = (1 / this.N) * np.fft.fftshift(np.fft.fft(u0.^2))

            # массивы с решениями   
            u = np.zeros((this.N, nt))
            u_hat = np.zeros((this.N, nt),  np.complex)
            u_hat2 = np.zeros((this.N, nt), np.complex)

            # инициализация массивов
            u[:, 1] = u0
            u_hat[:, 1] = u0_hat
            u_hat2[:, 1] = u0_hat2

            # основной цикл решателя
            for j in range(1,nt-1)
                uhat_curr = u_hat[:, j]
                uhat_curr2 = u_hat2[:,j]
                if j == 1
                    uhat_prev2 = u_hat2[:,1]
                else
                    uhat_prev2 = u_hat2[:,j-1]
                end
                # схема Кранка-Николсона для линейной части и Адамса-Башфорта для нелинейной части
                # таймстепинг в пространстве коэффициентов Фурье
                u_hat[:, j+1] .= ((1 ./ (1 .- (this.dt / 2) * this.FL)) .* ((1 .+ (this.dt / 2) * this.FL) .* uhat_curr + (((3 / 2) * this.FN) .* (uhat_curr2) - ((1 / 2) * this.FN) .* (uhat_prev2)) * this.dt))
                # возврат в физическое пространство
                u[:,j+1] = np.real(this.N * np.fft.ifft(np.fft.ifftshift(u_hat[:,j+1])))
                # корректировка коэффициента
                u_hat[:,j+1] = (1 / this.N) * np.fft.fftshift(np.fft.fft(u[:,j+1]))
                # вычисление нелинейной части
                u_hat2[:,j+1] = (1 / this.N) * np.fft.fftshift(np.fft.fft(u[:,j+1].^2))
            end

            this.u = u
        end

        this.show = function()
            fig, ax = plt.subplots(figsize=(25,6))
            tt, xx = np.meshgrid(this.t, this.x)
            levels = np.arange(-3, 3, 0.01)
            cs = ax.contourf(tt, xx, this.u)
            fig.colorbar(cs)
            ax.set_xlabel("t")
            ax.set_ylabel("x")
            ax.set_title("Kuramoto-Sivashinsky")
        end

    return this
    end
end

sol_ = Solver_(1, 35, 100, 0, 300, 0.05) # инициализируем решатель
sol_.solver() # решаем
u = sol_.u    # сохраняем решение
sol_.show()   # выводим решение