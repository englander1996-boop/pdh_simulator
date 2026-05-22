r"""column2/35.hsc のオブジェクト名を列挙する診断スクリプト。

partial condenser の column2 では Topl (液) と Topv (気) の両方が存在するはず。
adapter が正しいストリームを取れているか確認するために、HSC 内の
MaterialStreams / EnergyStreams / Operations を全列挙する。
"""
import os
import sys
import time

sys.path.insert(0, os.path.normpath(os.path.join(os.path.dirname(__file__), '..')))

from units.vle.hysys.registry import HsysRegistry
from units.vle.hysys.session import HysysSession


def main():
    reg = HsysRegistry()
    hsc_path = reg.get_path("column2", 35)
    print(f"HSC: {hsc_path}\n")

    with HysysSession(hsc_path, visible=False) as sess:
        fs = sess.flowsheet

        print("--- MaterialStreams ---")
        ms = fs.MaterialStreams
        for i in range(ms.Count):
            stream = ms.Item(i)
            name = stream.Name
            try:
                T = stream.Temperature.Value
                P = stream.Pressure.Value
                F = stream.MolarFlow.Value
                T_str = f"{T:.2f}" if abs(T - (-32767)) > 1 else "EMPTY"
                P_str = f"{P:.2f}" if abs(P - (-32767)) > 1 else "EMPTY"
                F_str = f"{F:.6f}" if abs(F - (-32767)) > 1 else "EMPTY"
                print(f"  [{i}] {name}: T={T_str} P={P_str} F={F_str}")
            except Exception as e:
                print(f"  [{i}] {name}: error {e}")

        print("\n--- EnergyStreams ---")
        es = fs.EnergyStreams
        for i in range(es.Count):
            stream = es.Item(i)
            try:
                Q = stream.HeatFlow.Value
                Q_str = f"{Q:.2f}" if abs(Q - (-32767)) > 1 else "EMPTY"
                print(f"  [{i}] {stream.Name}: Q={Q_str}")
            except Exception as e:
                print(f"  [{i}] {stream.Name}: error {e}")

        print("\n--- Operations ---")
        ops = fs.Operations
        for i in range(ops.Count):
            op = ops.Item(i)
            try:
                print(f"  [{i}] {op.Name} ({op.TypeName})")
            except Exception:
                print(f"  [{i}] {op.Name} (no TypeName)")

        # Tower (T-100) の内部構造を覗く
        print("\n--- T-100.ColumnFlowsheet ---")
        try:
            col = fs.Operations.Item("T-100")
            cf = col.ColumnFlowsheet

            print("  ColumnFlowsheet.MaterialStreams:")
            cfms = cf.MaterialStreams
            for i in range(cfms.Count):
                print(f"    [{i}] {cfms.Item(i).Name}")

            print("  ColumnFlowsheet.FeedStreams:")
            try:
                cffs = cf.FeedStreams
                for i in range(cffs.Count):
                    print(f"    [{i}] {cffs.Item(i).Name}")
            except Exception as e:
                print(f"    error: {e}")

            print("  ColumnFlowsheet.Operations:")
            cfo = cf.Operations
            for i in range(cfo.Count):
                op = cfo.Item(i)
                print(f"    [{i}] {op.Name}")

            print("  ColumnFlowsheet.Specifications (.Item):")
            try:
                specs = cf.Specifications
                for i in range(specs.Count):
                    s = specs.Item(i)
                    try:
                        cv = s.CurrentValue
                        gv = s.GoalValue
                        ia = s.IsActive
                        cv_str = f"{cv:.4f}" if abs(cv - (-32767)) > 1 else "EMPTY"
                        gv_str = f"{gv:.4f}" if abs(gv - (-32767)) > 1 else "EMPTY"
                        print(f"    [{i}] {s.Name}: Active={ia}, Cur={cv_str}, Goal={gv_str}")
                    except Exception as e:
                        print(f"    [{i}] {s.Name}: error {e}")
            except Exception as e:
                print(f"    Specifications error: {e}")
        except Exception as e:
            print(f"  T-100 error: {e}")


if __name__ == "__main__":
    main()
