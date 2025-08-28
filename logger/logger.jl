module MyLogger

using Logging
using Dates

export init_logger

struct FileLogger <: AbstractLogger
    # カスタム FileLogger をトップレベルで定義
    io::IO
    level::LogLevel
end

# ログレベルの判定
Logging.min_enabled_level(logger::FileLogger) = logger.level

function Logging.shouldlog(
        logger::FileLogger,
        level, _module, group, id
    )
    return level ≥ logger.level
end

function Logging.handle_message(
    logger::FileLogger,
    level, message, _module, group, id, file, line;
    kwargs...
    )
    # メッセージのハンドリング
    ts  = Dates.format(Dates.now(), dateformat"yyyy-mm-dd HH:MM:SS")
    lvl = uppercase(string(level))
    print(logger.io, "[", ts, "] [", lvl, "] ", message, "\n")
    flush(logger.io)
end

function init_logger(
    path::Union{Nothing,String}=nothing;
    dir::Union{Nothing,String}=nothing,
    level::Union{LogLevel,String,Symbol}=Logging.Info
    )
    """
        init_logger(
            path::Union{Nothing,String}=nothing;
            dir::Union{Nothing,String}=nothing,
            level::Union{LogLevel,String,Symbol}=Logging.Info
        )
    
    - `path`: ファイル名（省略時はタイムスタンプで自動生成）
    - `dir`: 保存先ディレクトリ（省略時はカレントディレクトリ）
    - `level`: ログレベルを `LogLevel`、文字列、またはシンボルで指定可能。無効なレベルを指定すると `ArgumentError` を投げる。
      - 例: `Logging.Debug`, `"Debug"`, `:Debug`, `"warning"`, `:ERROR` など。
    """
    # レベルの正規化
    lvl = begin
        if level isa LogLevel
            level
        elseif level isa Symbol || level isa String
            name = uppercase(string(level))
            mapping = Dict(
                "DEBUG" => Logging.Debug,
                "INFO" => Logging.Info,
                "WARN" => Logging.Warn,
                "WARNING" => Logging.Warn,
                "ERROR" => Logging.Error,
                "FATAL" => Logging.Error
            )
            if haskey(mapping, name)
                mapping[name]
            else
                throw(ArgumentError("Invalid log level: $level"))
            end
        else
            throw(ArgumentError("Invalid type for level: $(typeof(level))"))
        end
    end

    # ファイル名の決定
    filename = isnothing(path) ? string(Dates.format(Dates.now(), dateformat"yyyymmdd-HHMMSS"), ".log") : path

    # ディレクトリの作成（指定があれば）
    if dir !== nothing
        mkpath(dir)
    end

    # フルパスを組み立て
    filepath = dir === nothing ? filename : joinpath(dir, filename)

    # ファイルを追記モードで開く
    io = open(filepath, "a")
    # プログラム終了時にファイルを閉じる
    atexit(() -> close(io))

    # グローバルロガーに設定
    flogger = FileLogger(io, lvl)
    global_logger(flogger)

    @info "Logger initialized: path=$(filepath), level=$(lvl)"
    return nothing
end

end