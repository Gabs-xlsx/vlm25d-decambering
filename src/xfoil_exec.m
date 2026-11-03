function [st, expirou] = xfoil_exec(arq_cmd, arq_log, limite_s)
%XFOIL_EXEC  Roda "xfoil < arq_cmd > arq_log" COM LIMITE DE TEMPO.
%   [st, expirou] = xfoil_exec(cmd, log)          % limite padrao 300 s
%   [st, expirou] = xfoil_exec(cmd, log, 90)

    if nargin < 3 || isempty(limite_s); limite_s = 300; end
    limite_s = min(max(limite_s, 5), 1800);
    expirou  = false;
    linha    = sprintf('xfoil < "%s" > "%s"', arq_cmd, arq_log);

    lancou = false;
    if ispc
        try
            p = System.Diagnostics.Process();
            p.StartInfo.FileName         = 'cmd.exe';
            p.StartInfo.Arguments        = ['/c ' linha];
            p.StartInfo.UseShellExecute  = false;
            p.StartInfo.CreateNoWindow   = true;
            p.StartInfo.WorkingDirectory = pwd;
            p.Start();
            lancou = true;
            pid = double(p.Id);
            if p.WaitForExit(int32(round(limite_s*1000)))
                st = double(p.ExitCode);
            else
                expirou = true;  st = -1;
                [~,~] = system(sprintf('taskkill /F /T /PID %d > nul 2>&1', pid));
                try; p.WaitForExit(int32(3000)); catch; end
                if ~p.HasExited
                    [~,~] = system('taskkill /F /IM xfoil.exe > nul 2>&1');
                end
            end
            try; p.Close(); catch; end
            return
        catch ME
            if lancou
                [~,~] = system('taskkill /F /IM xfoil.exe > nul 2>&1');
                st = -1;  expirou = true;
                warning('xfoil:exec','falha ao esperar o XFOIL (%s); processo morto.', ME.message);
                return
            end
        end
    end

    [st,~] = system(linha);
end
